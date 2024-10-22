/*
 * SPDX-FileCopyrightText: 2024 Philip Abbet <philip.abbet@gmail.com>
 *
 * SPDX-FileContributor: Philip Abbet <philip.abbet@gmail.com>
 *
 * SPDX-License-Identifier: BSD-3-Clause
*/

#include <iostream>
#include <fstream>
#include <vector>
#include <filesystem>
#include <algorithm>

#define STB_IMAGE_IMPLEMENTATION
#include <stb_image.h>


extern "C" {
    #include <astrometry/simplexy.h>
    #include <astrometry/image2xy.h>
    #include <astrometry/index.h>
    #include <astrometry/solver.h>
    #include <astrometry/permutedsort.h>
}


/*********************************** HELPER FUNCTIONS ***********************************/

// Callback used to determine if a match is valid (always return true).
anbool match_callback(MatchObj* mo, void* userdata)
{
    return TRUE;
}

//-----------------------------------------------------------------------------

// Callback used to allow the application to cancel the plate solving.
// Here the first time it is called, we ask to be called back in 10 seconds,
// and the second time we stop the processing by returning 0.
time_t timer_callback(void* userdata)
{
    static bool mustStop = false;

    time_t limit = mustStop ? 0 : 10;
    mustStop = !mustStop;

    return limit;
}

//-----------------------------------------------------------------------------

// Load all the index files found in the given folder.
// Only the metadata of the files are really loaded at this stage.
std::vector<index_t*> loadIndexes(const std::string& folder)
{
    // Find all index files in the folder
    std::vector<std::string> files;

    if (!std::filesystem::is_directory(folder))
        return std::vector<index_t*>();

    for (const auto& dirEntry : std::filesystem::directory_iterator(folder))
    {
        if (dirEntry.is_directory())
            continue;

        if (dirEntry.path().extension() == ".fits")
        {
            if (index_is_file_index(dirEntry.path().c_str()))
                files.push_back(dirEntry.path().c_str());
        }
    }

    if (files.empty())
        return std::vector<index_t*>();

    // Sort them
    std::sort(files.begin(), files.end());

    // Load them in reverse
    std::vector<index_t*> indexes;
    for (auto iter = files.rbegin(); iter != files.rend(); ++iter)
    {
        index_t* index = index_load(iter->c_str(), INDEX_ONLY_LOAD_METADATA, nullptr);
        if (index)
            indexes.push_back(index);
    }

    return indexes;
}

//-----------------------------------------------------------------------------

// Returns the subset of the index files in 'index' that contain stars meaningful
// for the given 'minWidth' and 'maxWidth' (in degrees) and image size.
std::vector<index_t*> filterIndexes(
    const std::vector<index_t*>& indexes, float minWidth, float maxWidth,
    size_t imageWidth, size_t imageHeight
)
{
    if ((imageWidth == 0) || (imageHeight== 0))
        return std::vector<index_t*>();

    double scaleMin = deg2arcsec(minWidth) / imageWidth;
    double scaleMax = deg2arcsec(maxWidth) / imageWidth;

    // range of quad sizes that could be found in the field, in arcsec.
    double quadsize_min = 0.1 * std::min(imageWidth, imageHeight);
    double fmax = hypot(imageWidth, imageHeight) * scaleMax;
    double fmin = quadsize_min * scaleMin;

    std::vector<index_t*> result;
    for (index_t* index : indexes)
    {
        if (index_overlaps_scale_range(index, fmin, fmax))
            result.push_back(index);
    }

    return result;
}

//-----------------------------------------------------------------------------

// Returns a sorted list of star indices, that allows to retrieve star positions
// in the needed order.
std::vector<int> sort(const simplexy_t& params, bool ascending)
{
    if (params.npeaks <= 1)
        return std::vector<int>();

    int (*compare)(const void*, const void*);

    if (ascending)
        compare = compare_floats_asc;
    else
        compare = compare_floats_desc;

    // Set background = flux + background (ie, non-background-subtracted flux)
    float* background = new float[params.npeaks];
    for (int i = 0; i < params.npeaks; ++i)
        background[i] = params.background[i] + params.flux[i];

    // Sort by flux
    int* perm1 = permuted_sort(params.flux, sizeof(float), compare, NULL, params.npeaks);

    // Sort by non-background-subtracted flux
    int* perm2 = permuted_sort(background, sizeof(float), compare, NULL, params.npeaks);

    // Determine the final order
    bool* used = new bool[params.npeaks];
    memset(used, 0, params.npeaks * sizeof(bool));

    std::vector<int> result(params.npeaks);
    int k = 0;
    for (int i = 0; i < params.npeaks; ++i)
    {
        int inds[] = { perm1[i], perm2[i] };
        for (int j = 0; j < 2; ++j)
        {
            int index = inds[j];
            if (used[index])
                continue;

            used[index] = true;
            result[k] = index;
            ++k;
        }
    }

    delete[] background;
    delete[] used;
    free(perm1);
    free(perm2);

    return result;
}


/************************************* MAIN FUNCTION ************************************/

int main(int argc, char** argv)
{
    // Check the command-line arguments
    if (argc != 3)
    {
        std::cerr << "Usage: " << argv[0] << " <PNM image> <indexes folder>" << std::endl;
        return 1;
    }


    // Load the image as grayscale
    std::cout << "Loading the image '" << argv[1] << "'..." << std::endl;

    int imageWidth, imageHeight, imageChannels;
    unsigned char* pixels = stbi_load(argv[1], &imageWidth, &imageHeight, &imageChannels, 1);

    // Convert the image to floating-point
    std::vector<float> image(imageWidth * imageHeight);

    unsigned char* src = pixels;
    float* dst = image.data();

    for (int i = 0; i < imageWidth * imageHeight; ++i)
    {
        *dst = float(*src);
        ++src;
        ++dst;
    }

    stbi_image_free(pixels);


    // Detect stars
    std::cout << "Detecting stars..." << std::endl;

    simplexy_t params = { 0 };
    simplexy_fill_in_defaults(&params);

    params.image = image.data();
    params.nx = imageWidth;
    params.ny = imageHeight;

    int res = image2xy_run(&params, 2, 3);
    if (res != 0)
    {
        std::cerr << "Failed to detect stars" << std::endl;
        return 1;
    }

    std::cout << "   " << params.npeaks << " stars found" << std::endl;

    std::vector<int> sortedIndices = sort(params, false);
    
    params.image = nullptr;


    // Plate solving
    double minWidth = 0.1;
    double maxWidth = 180.0;
    // double minWidth = 0.5;
    // double maxWidth = 2.0;

    //--- load the index files
    std::cout << "Loading the index files from '" << argv[2] << "'..." << std::endl;

    std::vector<index_t*> indexes = loadIndexes(argv[2]);

    //--- create the solver
    solver_t* solver = solver_new();

    solver->pixel_xscale = 0.0;

    solver->field_maxx = double(imageWidth);
    solver->field_maxy = double(imageHeight);

    solver->funits_lower = deg2arcsec(0.1) / imageWidth;
    solver->funits_upper = deg2arcsec(180.0) / imageWidth;

    solver->logratio_toprint = log(1e6);
    solver->logratio_tokeep = log(1e9);
    solver->logratio_totune = log(1e6);

    solver->record_match_callback = match_callback;
    solver->timer_callback = timer_callback;
    solver->userdata = nullptr;

    solver->distance_from_quad_bonus = TRUE;
    solver->verify_dedup = FALSE;

    solver->do_tweak = TRUE;
    solver->tweak_aborder = 2;
    solver->tweak_abporder = 2;

    solver->quadsize_min = 0.1 * std::min(imageWidth, imageHeight);

    //--- only keep the relevant index files
    std::vector<index_t*> filteredIndexes = filterIndexes(indexes, minWidth, maxWidth, imageWidth, imageHeight);

    for (index_t* index : filteredIndexes)
    {
        if (!index->codekd)
        {
            if (index_reload(index) != 0)
                continue;
        }

        solver_add_index(solver, index);
    }

    std::cout << "    " << filteredIndexes.size() << " index files used" << std::endl;


    //--- copy the star positions
    starxy_t* fieldxy = starxy_new(std::min(params.npeaks, 1000), false, false);

    for (int i = 0; i < fieldxy->N; ++i)
    {
        int index = sortedIndices[i];

        fieldxy->x[i] = params.x[index];
        fieldxy->y[i] = params.y[index];
    }

    solver_set_field(solver, fieldxy);

    //--- now we can do plate solving
    std::cout << "Plate solving..." << std::endl;

    solver_run(solver);

    bool success = solver_did_solve(solver);
    if (success)
    {
        double ra,dec;
        xyzarr2radecdeg(solver->best_match.center, &ra, &dec);

        std::cout << "    " << ra << "°, " << dec << "°" << std::endl;
        std::cout << "    Pixel size: " << solver->best_match.scale << " arcsec" << std::endl;
    }
    else
    {
        std::cerr << "Failed to do plate solving" << std::endl;
        return 1;
    }

    solver_free(solver);
    simplexy_free_contents(&params);

    return 0;
}
