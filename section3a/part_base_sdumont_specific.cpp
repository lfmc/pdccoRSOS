// STD Libs
#include<cmath>
#include<string>
#include<ctime>
#include<cstdio>
#include<vector>

// VItA Libs
#include<structures/domain/AbstractDomain.h>
#include<structures/domain/SimpleDomain.h>
#include<structures/domain/DomainNVR.h>
#include<structures/domain/PartiallyVascularizedDomain.h>
#include<structures/domain/IntersectionVascularizedDomain.h>
#include<structures/domain/StagedDomain.h>
#include<core/GeneratorData.h>
#include<core/StagedFRROTreeGenerator.h>
#include<constrains/ConstantConstraintFunction.h>
#include<constrains/ConstantPiecewiseConstraintFunction.h>
#include<structures/tree/SingleVesselCCOOTree.h>
#include<structures/vascularElements/AbstractVascularElement.h>
#include<io/task/CheckpointSavingTask.h>
#include<io/task/VisualizationSavingTask.h>
#include<io/VTKObjectTreeNodalWriter.h>
#include<structures/tree/AdimSproutingVolumetricCostEstimator.h>
#include<structures/tree/SproutingVolumetricCostEstimator.h>
#include<core/TreeMerger.h>

#include<io/StagedFRROTreeGeneratorLogger.h>

using namespace std;

string base_folder {"/scratch/hmlffr2/luis.cury/"};
// string base_folder {"/home/lfmc/HeMoLab/VascularizacaoRenal/"};
string input_folder {base_folder + "Input/"};
string output_folder {base_folder + "Output/"};
string steps_folder {base_folder + "Steps/"};
string input_vtk_folder {input_folder + "VTK/"};
string input_cco_folder {input_folder + "CCO/"};
string output_vtp_folder {output_folder + "VTP/"};
string output_cco_folder {output_folder + "CCO/"};
string output_log_folder {output_folder + "LOG/"};
string output_points_folder {output_folder + "Points/"};
string steps_vtp_folder {steps_folder + "VTP/"};
string steps_cco_folder {steps_folder + "CCO/"};

string vascularize_base(point x0, string input_vtk, long long int total_term,vector<AbstractConstraintFunction<double, int>*> gam,
  vector<AbstractConstraintFunction<double, int>*> eps_lim, vector<AbstractConstraintFunction<double, int>*> nu,
  vector<int> n_draw, vector<int> seed, vector<int> n_level_test, vector<int> n_terminal_trial, vector<double> d_lim_red_factor,
  vector<double> perfusion_area_factor, vector<double> close_neighborhood_factor, vector<double> mid_point_d_lim_factor, vector<int> n_bifurcation_test,
  vector<double> min_bif_angle, double &last_Dlim, string suffix, double radius, double inflow)
{
    // Relevant input files
    //Geometry
    string input_vtk_hull {input_vtk_folder + input_vtk};
    // Output
    string output_base {"part_baseline_" + suffix + "base"};
    string output_cco {output_cco_folder + output_base + ".cco"};
    string output_log {output_log_folder + output_base + ".log"};
    string output_vtp {output_vtp_folder + output_base + ".vtp"};

    GeneratorData *gen_data_0 {new GeneratorData(n_level_test[0], n_terminal_trial[0], d_lim_red_factor[0], perfusion_area_factor[0], close_neighborhood_factor[0], mid_point_d_lim_factor[0], n_bifurcation_test[0], 0, false)};
    SimpleDomain *domain_0 {new SimpleDomain(input_vtk_hull, n_draw[0], seed[0], gen_data_0)};
    (*domain_0).setMinBifurcationAngle(min_bif_angle[0]);
    (*domain_0).setIsConvexDomain(true);

    GeneratorData *gen_data_1 {new GeneratorData(n_level_test[1], n_terminal_trial[1], d_lim_red_factor[1], perfusion_area_factor[1], close_neighborhood_factor[1], mid_point_d_lim_factor[1], n_bifurcation_test[1], 0, false)};
    SimpleDomain *domain_1 {new SimpleDomain(input_vtk_hull, n_draw[1], seed[1], gen_data_1)};
    (*domain_1).setMinBifurcationAngle(min_bif_angle[1]);
    (*domain_1).setIsConvexDomain(true);

    GeneratorData *gen_data_2 {new GeneratorData(n_level_test[2], n_terminal_trial[2], d_lim_red_factor[2], perfusion_area_factor[2], close_neighborhood_factor[2], mid_point_d_lim_factor[2], n_bifurcation_test[2], 0, false)};
    SimpleDomain *domain_2 {new SimpleDomain(input_vtk_hull, n_draw[2], seed[2], gen_data_2)};
    (*domain_2).setMinBifurcationAngle(min_bif_angle[2]);
    (*domain_2).setIsConvexDomain(true);

    lldiv_t term_division = lldiv(total_term, 6);
    vector<long long int> n_term {term_division.quot, 2 * term_division.quot, (3 * term_division.quot) + term_division.rem};

    StagedDomain *staged_domain {new StagedDomain()};
    (*staged_domain).addStage(n_term[0], domain_0);
    (*staged_domain).addStage(n_term[1], domain_1);
    (*staged_domain).addStage(n_term[2], domain_2);

    StagedFRROTreeGenerator *tree_generator {new StagedFRROTreeGenerator(staged_domain, x0, radius, inflow, total_term, gam, eps_lim, nu, 0.0, 1.e-6)};
    SingleVesselCCOOTree *tree  = static_cast<SingleVesselCCOOTree *> ((*tree_generator).getTree());
    (*tree).setIsInCm(true);

    (*tree_generator).generate(500, steps_folder);
    
    (*tree).save(output_cco);
    
    VTKObjectTreeNodalWriter *tree_writer {new VTKObjectTreeNodalWriter()};
    (*tree_writer).write(output_vtp, tree);
    
    FILE* fp {fopen(output_log.c_str(), "w")};
    if(!fp) {
        fprintf(stderr, "Failed to create configuration log file.\n");
        exit(EXIT_FAILURE);
    }
    StagedFRROTreeGeneratorLogger *logger {new StagedFRROTreeGeneratorLogger(fp, tree_generator)};
    (*logger).write();
    delete logger;
    fclose(fp);

    last_Dlim = (*tree_generator).getDLimLast();

    delete tree_writer;
    delete tree_generator;
    delete staged_domain;
    delete domain_2;
    delete gen_data_2;
    delete domain_1;
    delete gen_data_1;
    delete domain_0;
    delete gen_data_0;

    return output_cco;
}

string vascularize_part(string suffix, string input_part, string input_cco, long long int base_term, long long int total_term, double base_d_lim,
  vector<AbstractConstraintFunction<double, int>*> gam, vector<AbstractConstraintFunction<double, int>*> eps_lim, vector<AbstractConstraintFunction<double, int>*> nu,
  vector<int> n_draw, vector<int> seed, vector<int> n_level_test, vector<int> n_terminal_trial, vector<double> d_lim_red_factor,
  vector<double> perfusion_area_factor, vector<double> close_neighborhood_factor, vector<double> mid_point_d_lim_factor, vector<int> n_bifurcation_test,
  vector<double> min_bif_angle)
{
    // Relevant input files
    //Geometry
    string input_vtk_part {input_vtk_folder + input_part};
    // Output
    string output_base {"part_baseline_" + suffix};
    string output_cco {output_cco_folder + output_base + ".cco"};
    string output_log {output_log_folder + output_base + ".log"};
    string output_vtp {output_vtp_folder + output_base + ".vtp"};
    string output_points {output_points_folder + output_base + ".points"};

    GeneratorData *gen_data_0 {new GeneratorData(n_level_test[0], n_terminal_trial[0], d_lim_red_factor[0], perfusion_area_factor[0], close_neighborhood_factor[0], mid_point_d_lim_factor[0], n_bifurcation_test[0], 0, false)};
    SimpleDomain *domain_0 {new SimpleDomain(input_vtk_part, n_draw[0], seed[0], gen_data_0)};
    (*domain_0).setMinBifurcationAngle(min_bif_angle[0]);
    (*domain_0).setIsConvexDomain(false);

    GeneratorData *gen_data_1 {new GeneratorData(n_level_test[1], n_terminal_trial[1], d_lim_red_factor[1], perfusion_area_factor[1], close_neighborhood_factor[1], mid_point_d_lim_factor[1], n_bifurcation_test[1], 0, false)};
    SimpleDomain *domain_1 {new SimpleDomain(input_vtk_part, n_draw[1], seed[1], gen_data_1)};
    (*domain_1).setMinBifurcationAngle(min_bif_angle[1]);
    (*domain_1).setIsConvexDomain(false);

    lldiv_t term_division = lldiv(total_term, 3);
    vector<long long int> n_term {term_division.quot, 2 * term_division.quot + term_division.rem};

    StagedDomain *staged_domain {new StagedDomain()};
    (*staged_domain).addStage(n_term[0], domain_0);
    (*staged_domain).addStage(n_term[1], domain_1);
    int merge_stage {3};
    (*staged_domain).setInitialStage(merge_stage);

    SingleVesselCCOOTree *tree {new SingleVesselCCOOTree(input_cco, gen_data_0, gam[0], eps_lim[0], nu[0])};
    (*tree).setIsInCm(true);
    StagedFRROTreeGenerator *tree_generator {new StagedFRROTreeGenerator(staged_domain, tree, base_term + total_term, gam, eps_lim, nu)};
    (*tree_generator).setDLim(base_d_lim);

    FILE *fp_points {fopen(output_points.c_str(), "wb")};
    if(!fp_points) {
        fprintf(stderr, "Failed to create points file.\n");
        exit(EXIT_FAILURE);
    }
    (*tree_generator).resumeSavePointsMidPoint(10, steps_folder, fp_points);
    fclose(fp_points);

    
    (*tree).save(output_cco);
    VTKObjectTreeNodalWriter *tree_writer {new VTKObjectTreeNodalWriter()};
    (*tree_writer).write(output_vtp, tree);
    
    FILE* fp {fopen(output_log.c_str(), "w")};
    if(!fp) {
        fprintf(stderr, "Failed to create configuration log file.\n");
        exit(EXIT_FAILURE);
    }
    StagedFRROTreeGeneratorLogger *logger {new StagedFRROTreeGeneratorLogger(fp, tree_generator)};
    (*logger).write();
    delete logger;
    fclose(fp);
    
    delete tree_writer;
    delete tree_generator;
    delete tree;
    delete staged_domain;
    delete domain_1;
    delete gen_data_1;
    delete domain_0;
    delete gen_data_0;

    return output_points;
}

string vascularize_full_resume(string input_vtk, string input_cco, long long int base_term,long long int total_term,vector<AbstractConstraintFunction<double, int>*> gam,
  vector<AbstractConstraintFunction<double, int>*> eps_lim, vector<AbstractConstraintFunction<double, int>*> nu,
  vector<int> n_draw, vector<int> seed, vector<int> n_level_test, vector<int> n_terminal_trial, vector<double> d_lim_red_factor,
  vector<double> perfusion_area_factor, vector<double> close_neighborhood_factor, vector<double> mid_point_d_lim_factor, vector<int> n_bifurcation_test,
  vector<double> min_bif_angle, double &initial_Dlim, double &last_Dlim, string suffix)
{
    // Relevant input files
    //Geometry
    string input_vtk_hull {input_vtk_folder + input_vtk};
    // Output
    string output_base {"part_baseline_" + suffix  + "full"};
    string output_cco {output_cco_folder + output_base + ".cco"};
    string output_log {output_log_folder + output_base + ".log"};
    string output_vtp {output_vtp_folder + output_base + ".vtp"};

    GeneratorData *gen_data_0 {new GeneratorData(n_level_test[0], n_terminal_trial[0], d_lim_red_factor[0], perfusion_area_factor[0], close_neighborhood_factor[0], mid_point_d_lim_factor[0], n_bifurcation_test[0], 0, false)};
    SimpleDomain *domain_0 {new SimpleDomain(input_vtk_hull, n_draw[0], seed[0], gen_data_0)};
    (*domain_0).setMinBifurcationAngle(min_bif_angle[0]);
    (*domain_0).setIsConvexDomain(true);

    GeneratorData *gen_data_1 {new GeneratorData(n_level_test[1], n_terminal_trial[1], d_lim_red_factor[1], perfusion_area_factor[1], close_neighborhood_factor[1], mid_point_d_lim_factor[1], n_bifurcation_test[1], 0, false)};
    SimpleDomain *domain_1 {new SimpleDomain(input_vtk_hull, n_draw[1], seed[1], gen_data_1)};
    (*domain_0).setMinBifurcationAngle(min_bif_angle[1]);
    (*domain_0).setIsConvexDomain(true);

    lldiv_t term_division = lldiv(total_term, 3);
    vector<long long int> n_term {term_division.quot, 2 * term_division.quot + term_division.rem};

    StagedDomain *staged_domain {new StagedDomain()};
    (*staged_domain).addStage(n_term[0], domain_0);
    (*staged_domain).addStage(n_term[1], domain_1);
    (*staged_domain).setInitialStage(3);
    
    SingleVesselCCOOTree *tree = {new SingleVesselCCOOTree(input_cco, gen_data_0, gam[0], eps_lim[0], nu[0])};
    (*tree).setIsInCm(true);
    StagedFRROTreeGenerator *tree_generator {new StagedFRROTreeGenerator(staged_domain, tree, base_term + total_term, gam, eps_lim, nu)};
    (*tree_generator).setDLim(initial_Dlim);

    (*tree_generator).resume(500, steps_folder);
    
    (*tree).save(output_cco);
    
    VTKObjectTreeNodalWriter *tree_writer {new VTKObjectTreeNodalWriter()};
    (*tree_writer).write(output_vtp, tree);
    
    FILE* fp {fopen(output_log.c_str(), "w")};
    if(!fp) {
        fprintf(stderr, "Failed to create configuration log file.\n");
        exit(EXIT_FAILURE);
    }
    StagedFRROTreeGeneratorLogger *logger {new StagedFRROTreeGeneratorLogger(fp, tree_generator)};
    (*logger).write();
    delete logger;
    fclose(fp);

    last_Dlim = (*tree_generator).getDLimLast();

    delete tree_writer;
    delete tree_generator;
    delete tree;
    delete staged_domain;
    delete domain_0;
    delete gen_data_0;
    delete domain_1;
    delete gen_data_1;

    return output_cco;
}

void merge(string base_tree, vector<string> derived_points, AbstractConstraintFunction<double,int> *gam, AbstractConstraintFunction<double,int> *eps_lim,
    AbstractConstraintFunction<double,int> *nu, string suffix) {
    
    string output_base {"part_baseline_" + suffix + "merged"};
    string output_cco {output_cco_folder + output_base + ".cco"};
    string output_log {output_log_folder + output_base + ".log"};
    string output_vtp {output_vtp_folder + output_base + ".vtp"};

    GeneratorData *gen_data {new GeneratorData()};    
    SingleVesselCCOOTree *tree {new SingleVesselCCOOTree(base_tree, gen_data, gam, eps_lim, nu)};
    TreeMerger *tree_merger {new TreeMerger(tree, derived_points)};
    
    try {
        tree = {(*tree_merger).mergeFast()};
        (*tree).save(output_cco);
        VTKObjectTreeNodalWriter *tree_writer {new VTKObjectTreeNodalWriter()};
        (*tree_writer).write(output_vtp, tree);

        delete tree_writer;
    }
    catch(std::out_of_range const&) {
        std::cout << "This merge failed!" << std::endl;        
    }

    delete tree_merger;
    delete tree;
    delete gen_data;
}

int main(int argc, char *argv[])
{
    // if(argc != 2) {
    //     fprintf(stderr, "Did not include seed file.\n");
    //     exit(EXIT_FAILURE);
    // }

    string str = argv[1];
    char suffix_c[200];
    int base_seed_1, base_seed_2, base_seed_3;
    int region_seed_1, region_seed_2;
    FILE *fp = fopen(str.c_str(), "r");
    if(!fp) {
        fprintf(stderr, "Failed to open seed file.\n");
        exit(EXIT_FAILURE);
    }
    fscanf(fp, "%*s");
    fscanf(fp, "%d %d %d", &base_seed_1, &base_seed_2, &base_seed_3);
    fscanf(fp, "%d %d", &region_seed_1, &region_seed_2);
    fscanf(fp, "%s", suffix_c);
    vector<int> seed_base({base_seed_1, base_seed_2, base_seed_3});
    vector<int> seed_part({region_seed_1, region_seed_2});
    string suffix_base(suffix_c);
    int i {0};
    sscanf(argv[2], "%d", &i);
    fclose(fp);

    size_t max_no_stages {5};
    // double q0 {14.};
    int n_level_test {16000};
    int n_terminal_trial {2000};
    double d_lim_red_factor {.9};
    double mid_point_d_lim_factor {.25};
    int n_bif_test {7};
    AbstractConstraintFunction<double,int> *gam {new ConstantConstraintFunction<double, int>(3.)};
    AbstractConstraintFunction<double,int> *eps_lim {new ConstantPiecewiseConstraintFunction<double, int>({0.4, 0}, {0, 5})};
    // Viscosity in P
    AbstractConstraintFunction<double,int> *nu {new ConstantConstraintFunction<double, int>(3.6)}; //cP
    vector<int> n_draw_v(max_no_stages, 10000);
    vector<int> n_level_test_v(max_no_stages, n_level_test);
    vector<int> n_terminal_trial_v(max_no_stages, n_terminal_trial);
    vector<double> d_lim_red_factor_v(max_no_stages, d_lim_red_factor);
    vector<double> mid_point_d_lim_factor_v(max_no_stages, mid_point_d_lim_factor);
    vector<int> n_bif_test_v(max_no_stages, n_bif_test);
    vector<AbstractConstraintFunction<double,int> *> gam_v(max_no_stages, gam);
    vector<AbstractConstraintFunction<double,int> *> nu_v(max_no_stages, nu);
    vector<double> perfusion_area_factor_v(max_no_stages, 1.0);

    vector<AbstractConstraintFunction<double,int> *> eps_lim_v(max_no_stages, eps_lim);
    vector<double> min_bif_angle_v_full(max_no_stages, (3. / 18.) * M_PI);
    vector<double> close_neighbourhood_factor_v_full {8.0, 4.0, 2.0};
    
    vector<double> close_neighbourhood_factor_v_part {1.0, 0.5};

    string domain_full {"part_ar_d100.vtk"};

    vector<string> parts {"part_ar_d100_p01.vtk", "part_ar_d100_p02.vtk", "part_ar_d100_p03.vtk", "part_ar_d100_p04.vtk"};
    vector<string> suffix {"term0500_", "term1000_", "term1500_","term2000_", "term2500_", "term3000_", "term3500_", "term4000_", "term4500_"};

    vector<long long int> base_term {500, 1000, 1500, 2000, 2500, 3000, 3500, 4000, 4500};
    // vector<long long int> base_term(9, 100);
    vector<long long int> part_term {1125, 1000, 875, 750, 625, 500, 375, 250, 125};
    // vector<long long int> part_term(9, 500);

    double radius {0.15};
    double influx {7.};
    point input {0, -0.499, 0.1};

    
    double dlim_base {0};
    string base_cco = vascularize_base(input, domain_full, base_term[i], gam_v, eps_lim_v, nu_v, n_draw_v, seed_base, n_level_test_v, n_terminal_trial_v, d_lim_red_factor_v,
        perfusion_area_factor_v, close_neighbourhood_factor_v_full, mid_point_d_lim_factor_v, n_bif_test_v, min_bif_angle_v_full, dlim_base, suffix_base + "_" + suffix[i],
        radius, influx);
    vector<string> points;
    for (size_t j = 0; j < parts.size(); ++j) {
        char qualifier[30];
        sprintf(qualifier, "p%02lu_", j + 1);
        string qual_suffix = suffix_base + "_" + suffix[i] + string(qualifier);
        points.push_back(vascularize_part(qual_suffix, parts[j], base_cco, base_term[i], part_term[i], dlim_base, gam_v, eps_lim_v, nu_v,
            n_draw_v, seed_part, n_level_test_v, n_terminal_trial_v, d_lim_red_factor_v, perfusion_area_factor_v, close_neighbourhood_factor_v_part, mid_point_d_lim_factor_v,
            n_bif_test_v, min_bif_angle_v_full));
        }
    merge(base_cco, points, gam, eps_lim, nu, suffix_base + "_" + suffix[i]);
    double dlim_last {0};
    vascularize_full_resume(domain_full, base_cco, base_term[i], 4 * part_term[i], gam_v, eps_lim_v, nu_v, n_draw_v, seed_part, n_level_test_v,
        n_terminal_trial_v, d_lim_red_factor_v, perfusion_area_factor_v, close_neighbourhood_factor_v_part, mid_point_d_lim_factor_v, n_bif_test_v,
        min_bif_angle_v_full, dlim_base, dlim_last, suffix_base + "_" + suffix[i]);
    

    
    delete eps_lim;
    delete nu;
    delete gam;

    return 0;
} 
