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
#include<structures/tree/VolumetricCostEstimator.h>
#include<structures/tree/DoublePowerCostEstimator.h>
#include<structures/tree/LinearCombinationCostEstimator.h>
#include<structures/tree/PowerCostEstimator.h>
#include<core/TreeMerger.h>

#include<io/StagedFRROTreeGeneratorLogger.h>

using namespace std;

string base_folder {"/scratch/hmlffr2/luis.cury/"};
// string base_folder {"/home/lfmc/HeMoLab/VascularizacaoRenal/"};

string input_base_folder {base_folder + "Input/"};
string input_cco_folder {input_base_folder + "CCO/"};
string input_vtk_folder {input_base_folder + "VTK/"};

string output_base_folder {base_folder + "Output/"};
string output_cco_folder {output_base_folder + "CCO/"};
string output_vtp_folder {output_base_folder + "VTP/"};
string output_log_folder {output_base_folder + "LOG/"};
string output_points_folder {output_base_folder + "Points/"};

string steps_base_folder {base_folder + "Steps/"};
string steps_cco_folder {steps_base_folder + "CCO/"};
string steps_vtp_folder {steps_base_folder + "VTP/"};

string vascularize_cortex_region(string input_cco, int region, long long int region_terms, AbstractConstraintFunction<double, int> *gam, vector<int> seed, string suffix_id, double dlim, AbstractCostEstimator *cost)
{
    // Base geometry
    string input_vtk_surface {input_vtk_folder + "left_surface_06_07_2020_LFC_decimated.vtk"};
    string input_vtk_pelvis {input_vtk_folder + "left_pelvis_06_07_2020_LFC_decimated.vtk"};
    string input_vtk_pyramids {input_vtk_folder + "left_pyramids_06_07_2020_LFC_decimated.vtk"};
    string input_vtk_core {input_vtk_folder + "left_core_06_07_2020_LFC_decimated.vtk"};

    // This region
    char cortex_region_file[200];
    sprintf(cortex_region_file, "left_cortex_region_%02d_20_09_2020_cube.vtk", region);
    string input_vtk_cortex_region {input_vtk_folder + string(cortex_region_file)};

    char output_c_base[100];
    sprintf(output_c_base, "kidney_region_%02d_", region);
    string output_base {string(output_c_base) + suffix_id};

    int n_level_test {16000};
    int n_terminal_trial {2000};
    double d_lim_red_factor {.9};
    double mid_point_d_lim_factor {.25};
    double perfusion_area_factor {1.0};
    vector<double> close_neighborhood_factor({1.0, 0.5, 0.25, 0.125});
    int n_bifurcation_test {7};
    int n_draw {2000};
    vector<double> min_bif_angle(8, 3./18. * (M_PI));
    
    AbstractConstraintFunction<double, int> * eps_lim {new ConstantPiecewiseConstraintFunction<double, int>({0.4, 0}, {0, 7})};    
    AbstractConstraintFunction<double,int> *nu {new ConstantConstraintFunction<double, int>(3.6)}; //cP
    vector<AbstractConstraintFunction<double,int> *> gams(8, gam);
    vector<AbstractConstraintFunction<double,int> *> eps_lims(8, eps_lim);
    vector<AbstractConstraintFunction<double,int> *> nus(8, nu);

    long long base_terms {5000};
    lldiv_t term_division = lldiv(region_terms, 10);
    long long int n_term[4] {term_division.quot, 2 * term_division.quot, 3 * term_division.quot,  (4 * term_division.quot) + term_division.rem};

    GeneratorData *gen_data_0 {new GeneratorData(n_level_test, n_terminal_trial, d_lim_red_factor, perfusion_area_factor, close_neighborhood_factor[0], mid_point_d_lim_factor, n_bifurcation_test, 0, false, cost)};
    PartiallyVascularizedDomain *domain_0 {new PartiallyVascularizedDomain(input_vtk_cortex_region, {input_vtk_surface}, {input_vtk_core, input_vtk_pelvis, input_vtk_pyramids}, n_draw, seed[0], gen_data_0)};
    (*domain_0).setMinBifurcationAngle(min_bif_angle[0]);
    (*domain_0).setIsConvexDomain(false);

    GeneratorData *gen_data_1 {new GeneratorData(n_level_test, n_terminal_trial, d_lim_red_factor, perfusion_area_factor, close_neighborhood_factor[1], mid_point_d_lim_factor, n_bifurcation_test, 0, false, cost)};
    PartiallyVascularizedDomain *domain_1 {new PartiallyVascularizedDomain(input_vtk_cortex_region, {input_vtk_surface}, {input_vtk_core, input_vtk_pelvis, input_vtk_pyramids}, n_draw, seed[1], gen_data_1)};
    (*domain_1).setMinBifurcationAngle(min_bif_angle[1]);
    (*domain_1).setIsConvexDomain(false);

    GeneratorData *gen_data_2 {new GeneratorData(n_level_test, n_terminal_trial, d_lim_red_factor, perfusion_area_factor, close_neighborhood_factor[2], mid_point_d_lim_factor, n_bifurcation_test, 0, false, cost)};
    PartiallyVascularizedDomain *domain_2 {new PartiallyVascularizedDomain(input_vtk_cortex_region, {input_vtk_surface}, {input_vtk_core, input_vtk_pelvis, input_vtk_pyramids}, n_draw, seed[2], gen_data_2)};
    (*domain_2).setMinBifurcationAngle(min_bif_angle[2]);
    (*domain_2).setIsConvexDomain(false);

    GeneratorData *gen_data_3 {new GeneratorData(n_level_test, n_terminal_trial, d_lim_red_factor, perfusion_area_factor, close_neighborhood_factor[3], mid_point_d_lim_factor, n_bifurcation_test, 0, false, cost)};
    PartiallyVascularizedDomain *domain_3 {new PartiallyVascularizedDomain(input_vtk_cortex_region, {input_vtk_surface}, {input_vtk_core, input_vtk_pelvis, input_vtk_pyramids}, n_draw, seed[3], gen_data_3)};
    (*domain_3).setMinBifurcationAngle(min_bif_angle[3]);
    (*domain_3).setIsConvexDomain(false);

    StagedDomain *staged_domain {new StagedDomain()};
    (*staged_domain).addStage(n_term[0], domain_0);
    (*staged_domain).addStage(n_term[1], domain_1);
    (*staged_domain).addStage(n_term[2], domain_2);
    (*staged_domain).addStage(n_term[3], domain_3);
    int merge_stage {4};
    (*staged_domain).setInitialStage(merge_stage);

    SingleVesselCCOOTree *tree {new SingleVesselCCOOTree(input_cco, gen_data_0, gam, eps_lim, nu)};
    (*tree).setIsInCm(true);
    StagedFRROTreeGenerator *tree_generator {new StagedFRROTreeGenerator(staged_domain, tree, base_terms + region_terms, gams, eps_lims, nus)};
    (*tree_generator).setDLim(dlim);
    CheckpointSavingTask *cco_saving {new CheckpointSavingTask(steps_cco_folder, "_step")};
    (*tree_generator).setSavingTasks({cco_saving});
       
    string output_points = output_points_folder + output_base + ".points";
    FILE *fp_points = fopen(output_points.c_str(), "wb");
    if (!fp_points) {
        fprintf(stderr, "Failed to create points files!\n");
        exit(EXIT_FAILURE);
    }
    tree = {(SingleVesselCCOOTree *) (*tree_generator).resumeSavePointsMidPoint(300, steps_base_folder, fp_points)};
    fclose(fp_points);
    
    (*tree).save(output_cco_folder + output_base + ".cco");
    VTKObjectTreeNodalWriter *tree_writer {new VTKObjectTreeNodalWriter()};
    (*tree_writer).write(output_vtp_folder + output_base + ".vtp", tree);
    
    FILE* fp {fopen((output_log_folder + output_base + ".log").c_str(), "w")};
    if(!fp) {
        fprintf(stderr, "Failed to create configuration log file.\n");
        exit(EXIT_FAILURE);
    }
    StagedFRROTreeGeneratorLogger *logger {new StagedFRROTreeGeneratorLogger(fp, tree_generator)};
    (*logger).write();
    delete logger;
    fclose(fp);
      
    delete cco_saving;
    delete tree_writer;
    delete tree_generator;
    delete tree;
    delete staged_domain;
    delete domain_3;
    delete gen_data_3;
    delete domain_2;
    delete gen_data_2;
    delete domain_1;
    delete gen_data_1;
    delete domain_0;
    delete gen_data_0;

    delete eps_lim;
    delete nu;

    return output_points;
}

AbstractCostEstimator * read_cost(FILE *fp) {
    char cost_type[300];
    fscanf(fp, "%s\n", cost_type);
    if (!strcmp(cost_type, "VolumetricCostEstimator")) {
        return (new VolumetricCostEstimator());
    }
    else if (!strcmp(cost_type, "SproutingVolumetricCostEstimator")) {
        double vf, pf, df;
        fscanf(fp, "%lf %lf %lf\n", &vf, &pf, &df);
        return (new SproutingVolumetricCostEstimator(vf, pf, df));
    }
    else if (!strcmp(cost_type, "AdimSproutingVolumetricCostEstimator")) {
        double vf, pf, df, vref, rref;
        fscanf(fp, "%lf %lf %lf %lf %lf\n", &vf, &pf, &df, &vref, &rref);
        return (new AdimSproutingVolumetricCostEstimator(vf, pf, df, vref, rref));
    }
    else if (!strcmp(cost_type, "PowerCostEstimator")) {
        unsigned int l_exp, r_exp;
        fscanf(fp, "%u %u\n", &l_exp, &r_exp);
        return (new PowerCostEstimator(l_exp, r_exp));
    }
    else if (!strcmp(cost_type, "DoublePowerCostEstimator")) {
        double l_exp, r_exp;
        fscanf(fp, "%lf %lf\n", &l_exp, &r_exp);
        return (new DoublePowerCostEstimator(l_exp, r_exp));
    }
    else if (!strcmp(cost_type, "LinearCombinationCostEstimator")) {
        unsigned int size;
        vector<double> alphas;
        vector<unsigned int> l_exps, r_exps; 
        fscanf(fp, "%u\n", &size);
        for (unsigned int i = 0; i < size; ++i) {
            double alpha;
            unsigned int l_exp, r_exp;
            fscanf(fp, "%lf %u %u\n", &alpha, &l_exp, &r_exp);
            alphas.push_back(alpha);
            l_exps.push_back(l_exp);
            r_exps.push_back(r_exp);
        }
        return (new LinearCombinationCostEstimator(alphas, l_exps, r_exps));
    }
    else {
        fprintf(stderr, "Bad cost estimator configuration.\n");
        exit(EXIT_FAILURE);
    }
}

int main(int argc, char *argv[])
{   
    if (argc != 3) {
        fprintf(stderr, "exec conf_file part_no");
        exit(EXIT_FAILURE);
    }
    string str = argv[1];
    int part_no {0};
    sscanf(argv[2], "%d", &part_no);

    int part_seed_1, part_seed_2, part_seed_3, part_seed_4;
    char suffix_c[100];
    double gamma_value;

    double dlim;
    
    AbstractCostEstimator *cost_functional {nullptr};
    AbstractConstraintFunction<double, int> *gam {nullptr};

    FILE *fp = fopen(str.c_str(), "r");
    if(!fp) {
        fprintf(stderr, "Failed to open seed file.\n");
        exit(EXIT_FAILURE);
    }
    fscanf(fp, "%d %d %d %d\n", &part_seed_1, &part_seed_2, &part_seed_3, &part_seed_4);
    vector<int> seed {part_seed_1, part_seed_2, part_seed_3, part_seed_4};
    fscanf(fp, "%lf\n", &dlim);
    cost_functional = read_cost(fp);
    fscanf(fp, "%lf\n", &gamma_value);
    gam = new ConstantConstraintFunction<double, int>(gamma_value);
    fscanf(fp, "%s\n", suffix_c);
    printf("%s\n", suffix_c);
    fclose(fp);

    string input_cco {input_cco_folder + "kidney_base_" + string(suffix_c) + ".cco"};
    // printf("input_cco = %s\n", input_cco.c_str());
    
    // Terms for 50k terminals kidney
    // vector<long long int> term_per_region {3912, 7007, 5964, 6033, 6973, 6546, 7413, 4152};
    // Terms for 100k terminals kidney
    vector<long long int> term_per_region {7742, 13870 , 11809, 11942, 13803, 12958 , 14668, 8208};

    vascularize_cortex_region(input_cco, part_no, term_per_region[part_no - 1], gam, seed, string(suffix_c), dlim, cost_functional);

    delete gam;
    delete cost_functional;

    return 0;
} 
