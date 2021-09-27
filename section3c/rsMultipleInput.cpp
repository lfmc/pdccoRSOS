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
#include<core/TreeMerger.h>
#include<io/task/CheckpointSavingTask.h>
#include<io/task/VisualizationSavingTask.h>
#include<io/VTKObjectTreeNodalWriter.h>
#include<io/VTKObjectTreeStrahlerWriter.h>
#include<structures/tree/AdimSproutingVolumetricCostEstimator.h>
#include<structures/tree/SproutingVolumetricCostEstimator.h>

#include<io/StagedFRROTreeGeneratorLogger.h>

using namespace std;

// Folders location
string input_folder {"/home/lfmc/HeMoLab/VascularizacaoRenal/Input/"};
string output_folder {"/home/lfmc/HeMoLab/VascularizacaoRenal/Output/"};
string steps_folder {"/home/lfmc/HeMoLab/VascularizacaoRenal/Steps/"};
string input_vtk_folder {input_folder + "VTK/"};
string input_cco_folder {input_folder + "CCO/"};
string output_vtp_folder {output_folder + "VTP/"};
string output_cco_folder {output_folder + "CCO/"};
string output_log_folder {output_folder + "LOG/"};
string output_points_folder {output_folder + "Points/"};
string steps_vtp_folder {steps_folder + "VTP/"};
string steps_cco_folder {steps_folder + "CCO/"};

string vascularize_base(int vol_f, int prot_f, int dif_f, string suffix_id)
{
    printf("Begin base tree vascularization.\n");

    double v0 {8};
    double r0 {0.09};
    double q0 {0.72};
    int n_level_test {16000};
    int n_terminal_trial {2000};
    double d_lim_red_factor {.9};
    double mid_point_d_lim_factor {.25};
    double perfusion_area_factor {1.0};
    double close_neighborhood_factor {8.0};
    int n_bifurcation_test {7};
    int n_draw {2000};
    int seed {1922};
    double min_bif_angle {3./18. * (M_PI)};
    long long int n_term {200};
    
    AbstractConstraintFunction<double, int> *gam {new ConstantConstraintFunction<double, int>(3.0)};
    AbstractConstraintFunction<double, int> *nu {new ConstantConstraintFunction<double, int>(3.6)};
    AbstractConstraintFunction<double, int> *eps_lim {new ConstantPiecewiseConstraintFunction<double, int>({0.4, 0}, {0, 7})};

    // Timestamp   
    // time_t now_time_t = time(nullptr);
    // struct tm *now_tm = localtime(&now_time_t);
    // char time_c_string[21];
    // strftime(time_c_string, 20, "%d_%m_%Y_%H_%M_%S", now_tm);
    // string time_string(time_c_string);


    // Relevant input files
    //Geometry
    string input_vtk_hull {input_vtk_folder + "rsMultipleInputhull.vtk"};
    //Vascular tree
    string input_cco {input_cco_folder + "rsMultipleInput_new.cco"};
    // Output
    
    // string output_base {"rsMultipleInputBase_" + string(functional_name) + "_" + time_string};
    string output_base {"rsMultipleInput_" + suffix_id};
    
    double scale_factor {10000.0};
	AdimSproutingVolumetricCostEstimator *cost_estimator {new AdimSproutingVolumetricCostEstimator(vol_f / scale_factor, prot_f / scale_factor, dif_f / scale_factor, v0, r0)};
    AbstractVascularElement::VESSEL_FUNCTION function {AbstractVascularElement::VESSEL_FUNCTION::DISTRIBUTION};
    GeneratorData *gen_data_0 {new GeneratorData(n_level_test, n_terminal_trial, 
     d_lim_red_factor, perfusion_area_factor, close_neighborhood_factor, mid_point_d_lim_factor, 
     n_bifurcation_test, function, false, cost_estimator)};

    SimpleDomain *domain_0 {new SimpleDomain(input_vtk_hull, n_draw, seed, gen_data_0)};

    (*domain_0).setIsConvexDomain(false);
    (*domain_0).setMinBifurcationAngle(min_bif_angle);
    StagedDomain *staged_domain {new StagedDomain()};
    (*staged_domain).addStage(n_term, domain_0);

    long long int n_term_total {n_term};

    SingleVesselCCOOTree *tree {new SingleVesselCCOOTree(input_cco, gen_data_0, q0, gam, eps_lim, nu, 0.0, 1.0e-5)};
    tree->setIsInCm(true);
    StagedFRROTreeGenerator *tree_generator {new StagedFRROTreeGenerator(staged_domain, tree, n_term_total, {gam}, {eps_lim}, {nu})};

    // CheckpointSavingTask *check_saving_task {new CheckpointSavingTask(steps_folder, output_base + "_step_")};
    // VisualizationSavingTask *vtk_saving_task {new VisualizationSavingTask(steps_folder, output_base + "_steps_")};

    // (*tree_generator).setSavingTasks({check_saving_task, vtk_saving_task});

    tree = {(SingleVesselCCOOTree *) (*tree_generator).resume(500, steps_folder)};
    string output_cco {output_cco_folder + output_base + ".cco"};
    (*tree).save(output_cco);
    VTKObjectTreeStrahlerWriter *tree_writer {new VTKObjectTreeStrahlerWriter()};
    string output_vtp {output_vtp_folder + output_base + ".vtp"};
    (*tree_writer).write(output_vtp, tree);

    string output_log {output_log_folder + output_base + ".log"};
    FILE* fp {fopen((output_log).c_str(), "w")};
    if(!fp) {
        fprintf(stderr, "Failed to create configuration log file.\n");
        exit(EXIT_FAILURE);
    }
    StagedFRROTreeGeneratorLogger *logger {new StagedFRROTreeGeneratorLogger(fp, tree_generator)};
    (*logger).write();
    delete logger;
    fclose(fp);
    
    delete tree_writer;
    // delete vtk_saving_task;
    // delete check_saving_task;
    // delete tree_generator;
    delete tree;
    delete staged_domain;
    delete domain_0;
    delete gen_data_0;
    delete eps_lim;
    delete nu;
    delete gam;

    printf("End base tree vascularization");
    return output_cco;
}

string vascularize_part(string suffix, string input_part, string input_cco)
{
    // Relevant input files
    //Geometry
    string input_vtk_part {input_vtk_folder + input_part};
    // Output
    string output_base {"rsMultipleInput_" + suffix};
    string output_cco {output_cco_folder + output_base + ".cco"};
    string output_log {output_log_folder + output_base + ".log"};
    string output_vtp {output_vtp_folder + output_base + ".vtp"};
    string output_points {output_points_folder + output_base + ".points"};

    int n_level_test {16000};
    int n_terminal_trial {2000};
    double d_lim_red_factor {.9};
    double mid_point_d_lim_factor {.25};
    double perfusion_area_factor {1.0};
    vector<double> close_neighborhood_factor ({4., 2.});
    int n_bifurcation_test {7};
    int n_draw {2000};
    vector<int> seed ({1517, 32});
    double min_bif_angle {3./18. * (M_PI)};
    long long int base_term {200};
    long long int part_term {400};
    AbstractVascularElement::VESSEL_FUNCTION function {AbstractVascularElement::VESSEL_FUNCTION::DISTRIBUTION};
    
    AbstractConstraintFunction<double, int> *gam {new ConstantConstraintFunction<double, int>(3.0)};
    AbstractConstraintFunction<double, int> *nu {new ConstantConstraintFunction<double, int>(3.6)};
    AbstractConstraintFunction<double, int> *eps_lim {new ConstantPiecewiseConstraintFunction<double, int>({0.4, 0}, {0, 7})};



    GeneratorData *gen_data_0 {new GeneratorData(n_level_test, n_terminal_trial, d_lim_red_factor, perfusion_area_factor, close_neighborhood_factor[0], mid_point_d_lim_factor, n_bifurcation_test, function, false)};
    SimpleDomain *domain_0 {new SimpleDomain(input_vtk_part, n_draw, seed[0], gen_data_0)};
    (*domain_0).setMinBifurcationAngle(min_bif_angle);
    (*domain_0).setIsConvexDomain(false);

    GeneratorData *gen_data_1 {new GeneratorData(n_level_test, n_terminal_trial, d_lim_red_factor, perfusion_area_factor, close_neighborhood_factor[1], mid_point_d_lim_factor, n_bifurcation_test, function, false)};
    SimpleDomain *domain_1 {new SimpleDomain(input_vtk_part, n_draw, seed[1], gen_data_1)};
    (*domain_1).setMinBifurcationAngle(min_bif_angle);
    (*domain_1).setIsConvexDomain(false);

    lldiv_t term_division = lldiv(part_term, 3);
    vector<long long int> n_term {term_division.quot, 2 * term_division.quot + term_division.rem};

    StagedDomain *staged_domain {new StagedDomain()};
    (*staged_domain).addStage(n_term[0], domain_0);
    (*staged_domain).addStage(n_term[1], domain_1);
    int merge_stage {1};
    (*staged_domain).setInitialStage(merge_stage);

    SingleVesselCCOOTree *tree {new SingleVesselCCOOTree(input_cco, gen_data_0, gam, eps_lim, nu)};
    (*tree).setIsInCm(true);
    StagedFRROTreeGenerator *tree_generator {new StagedFRROTreeGenerator(staged_domain, tree, base_term + part_term, {gam, gam, gam}, {eps_lim, eps_lim, eps_lim}, {nu, nu, nu})};
    // (*tree_generator).setDLim(base_d_lim);

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

    delete eps_lim;
    delete nu;
    delete gam;

    return output_points;
}

void merge(string base_tree, vector<string> derived_points, string suffix) {
    
    printf("Begin merge.\n");
    string output_base {"rsMultipleInput_" + suffix };
    string output_cco {output_cco_folder + output_base + ".cco"};
    string output_log {output_log_folder + output_base + ".log"};
    string output_vtp {output_vtp_folder + output_base + ".vtp"};

    AbstractConstraintFunction<double, int> *gam {new ConstantConstraintFunction<double, int>(3.0)};
    AbstractConstraintFunction<double, int> *nu {new ConstantConstraintFunction<double, int>(3.6)};
    AbstractConstraintFunction<double, int> *eps_lim {new ConstantPiecewiseConstraintFunction<double, int>({0.4, 0}, {0, 7})};

    GeneratorData *gen_data {new GeneratorData()};
    SingleVesselCCOOTree *tree {new SingleVesselCCOOTree(base_tree, gen_data, gam, eps_lim, nu)};
    TreeMerger *tree_merger {new TreeMerger(tree, derived_points)};
    
    try {
        tree = {(*tree_merger).mergeFast()};
        (*tree).save(output_cco);
        VTKObjectTreeStrahlerWriter *tree_writer {new VTKObjectTreeStrahlerWriter()};
        (*tree_writer).write(output_vtp, tree);

        delete tree_writer;
    }
    catch(std::out_of_range const&) {
        std::cout << "This merge failed!" << std::endl; 
    }

    delete tree_merger;
    delete tree;
    delete gen_data;
    delete eps_lim;
    delete nu;
    delete gam;

    printf("End merge.\n");
}

int main(int argc, char *argv[])
{
    // int factors[4][3] {{9990, 5, 5}, {9990, 10, 0}, {9990, 0, 10}, {10000, 0, 0}};
    // int factors[10][3] {{2000, 0, 8000}, {1000, 0, 9000}, {3000, 0, 7000}, 
    //     {8000, 0, 2000}, {6000, 0, 4000}, {4000, 0, 6000}, {9000, 0, 1000}, {7000, 0, 3000},
    //     {5000, 0, 5000}, {0, 0, 10000}};
    int factors[3][3] {{1000, 8000, 1000}, {2000, 7000, 1000}, {4000, 5000, 1000}};
    vector<string> parts {"rsMultipleInputp01.vtk", "rsMultipleInputp02.vtk", "rsMultipleInputp03.vtk", 
        "rsMultipleInputp04.vtk", "rsMultipleInputp05.vtk", "rsMultipleInputp06.vtk", "rsMultipleInputp07.vtk", 
        "rsMultipleInputp08.vtk"};

    // for (size_t vol_idx = 0; vol_idx < 11; ++vol_idx) {
    //     for (size_t prot_idx = 0; prot_idx < vol_idx; ++prot_idx) {
    //         int vol_f = 10000 - (1000 * vol_idx);
    //         int prot_f = 1000 * prot_idx;
    //         int dif_f = 10000 - vol_f - prot_f;
    //         char functional_name[200];
    //         sprintf(functional_name, "func_%05d_%05d_%05d", vol_f, prot_f, dif_f);
    //         string cco_base {vascularize_base(vol_f, prot_f, dif_f, string(functional_name)+"_base")};
    //     }
    for (int i = 0 ; i < 3; ++i) {
        int vol_f, prot_f, dif_f;
        vol_f = factors[i][0];
        prot_f = factors[i][1];
        dif_f = factors[i][2];
        char functional_name[200];
        sprintf(functional_name, "func_%05d_%05d_%05d", vol_f, prot_f, dif_f);
        // string cco_base {vascularize_base(vol_f, prot_f, dif_f, string(functional_name)+"_base")};
        string cco_base {output_cco_folder + "rsMultipleInput_" + string(functional_name)+"_base.cco"};
        vector<string> points;
        for (size_t j = 0; j != parts.size(); ++j) {
            char part_name[50];
            sprintf(part_name, "_p%02lu", j);
            string points_output {vascularize_part(string(functional_name)+ string(part_name), parts[j], cco_base)};
            points.push_back(points_output);
        }
        merge(cco_base, points, string(functional_name) + "_merged");
    }    

    return 0;
}