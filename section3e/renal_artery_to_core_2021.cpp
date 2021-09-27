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

#include<io/StagedFRROTreeGeneratorLogger.h>

using namespace std;

void set_distal_terminal(SingleVessel *root) 
{
    if (!root) {
        return;
    }
    vector<AbstractVascularElement *> children = (*root).getChildren();
    if(children.empty()) {
        (*root).branchingMode = AbstractVascularElement::BRANCHING_MODE::DISTAL_BRANCHING;
    }
    else {
        (*root).branchingMode = AbstractVascularElement::BRANCHING_MODE::NO_BRANCHING;
        for(auto vessel = children.begin(); vessel != children.end(); ++vessel) {
            set_distal_terminal((SingleVessel *) *vessel);
        }
    }

}

void make_terminal(SingleVesselCCOOTree *tree) 
{
    set_distal_terminal((SingleVessel *) (*tree).getRoot());
}

void adan_kidney_vascularize_cco_functional(double q0, vector<long long int> n_term,vector<AbstractConstraintFunction<double, int>*> gam,
  vector<AbstractConstraintFunction<double, int>*> eps_lim, vector<AbstractConstraintFunction<double, int>*> nu,
  vector<int> n_draw, vector<int> seed, vector<int> n_level_test, vector<int> n_terminal_trial, vector<double> d_lim_red_factor,
  vector<double> perfusion_area_factor, vector<double> close_neighborhood_factor, vector<double> mid_point_d_lim_factor, vector<int> n_bifurcation_test,
  vector<double> min_bif_angle, int vol_ref, int rad_ref, int vol_f, int prot_f, int dif_f, string cco_input, string suffix_id)
{
    printf("Beginning adan_kidney_vascularize.\n");

    // Timestamp   
    time_t now_time_t = time(nullptr);
    struct tm *now_tm = localtime(&now_time_t);
    char time_c_string[21];
    strftime(time_c_string, 20, "%d_%m_%Y_%H_%M_%S", now_tm);
    string time_string(time_c_string);

    // Folders location
    string input_folder {"/home/lfmc/HeMoLab/VascularizacaoRenal/Input/"};
    string output_folder {"/home/lfmc/HeMoLab/VascularizacaoRenal/Output/"};
    string steps_folder {"/home/lfmc/HeMoLab/VascularizacaoRenal/Steps/"};
    // Relevant input files
    //Geometry
    string input_vtk_surface {input_folder + "VTK/left_surface_06_07_2020_LFC_decimated.vtk"};
    string input_vtk_pelvis {input_folder + "VTK/left_pelvis_06_07_2020_LFC_decimated.vtk"};
    string input_vtk_pyramids {input_folder + "VTK/left_pyramids_06_07_2020_LFC_decimated.vtk"};
    string input_vtk_core {input_folder + "VTK/left_core_06_07_2020_LFC_decimated.vtk"};
    string input_vtk_cortex {input_folder + "VTK/left_cortex_06_07_2020_LFC_decimated.vtk"};
    string input_vtk_cube {input_folder + "VTK/test_cube_v2_09_07_2020.vtk"};
    //Vascular tree
    string input_cco {input_folder + "CCO/" + cco_input};
    // Output
    char functional_name[200];
    sprintf(functional_name, "func_%05d_%05d_%05d", vol_f, prot_f, dif_f);
    string output_base {"left_kidney_13_pyramids_" + suffix_id + "_" + string(functional_name) + "_" + time_string};
    
    
    double scale_factor {10000.0};
	AdimSproutingVolumetricCostEstimator *cost_estimator {new AdimSproutingVolumetricCostEstimator(vol_f / scale_factor, prot_f / scale_factor, dif_f / scale_factor, vol_ref, rad_ref)};
    AbstractVascularElement::VESSEL_FUNCTION function {AbstractVascularElement::VESSEL_FUNCTION::PERFORATOR};
    GeneratorData *gen_data_0 {new GeneratorData(n_level_test[0], n_terminal_trial[0], d_lim_red_factor[0], perfusion_area_factor[0], close_neighborhood_factor[0], mid_point_d_lim_factor[0], n_bifurcation_test[0], function, false, cost_estimator)};
    DomainNVR *domain_0 {new DomainNVR(input_vtk_core, {input_vtk_pelvis, input_vtk_pyramids}, n_draw[0], seed[0], gen_data_0)};
    (*domain_0).setIsConvexDomain(false);
    (*domain_0).setMinBifurcationAngle(min_bif_angle[0]);


    StagedDomain *staged_domain {new StagedDomain()};
    (*staged_domain).addStage(n_term[0], domain_0);

    long long int n_term_total {n_term[0]};

    SingleVesselCCOOTree *tree {new SingleVesselCCOOTree(input_cco, gen_data_0, q0, gam[0], eps_lim[0], nu[0], 0.0, 1.0e-6)};
    tree->setIsInCm(true);
    StagedFRROTreeGenerator *tree_generator {new StagedFRROTreeGenerator(staged_domain, tree, n_term_total, gam, eps_lim, nu)};

    tree = {(SingleVesselCCOOTree *) (*tree_generator).resume(500, steps_folder)};
    make_terminal(tree);
    if(!tree) {
        printf("Failed to generate tree!\n");
        return;
    }
    (*tree).save(output_folder + "CCO/" + output_base + ".cco");
    VTKObjectTreeNodalWriter *tree_writer {new VTKObjectTreeNodalWriter()};
    (*tree_writer).write(output_folder + "VTP/" + output_base + ".vtp", tree);

    FILE* fp {fopen((output_folder + "LOG/" + output_base + ".log").c_str(), "w")};
    if(!fp) {
        fprintf(stderr, "Failed to create configuration log file.\n");
        exit(EXIT_FAILURE);
    }
    StagedFRROTreeGeneratorLogger *logger {new StagedFRROTreeGeneratorLogger(fp, tree_generator)};
    (*logger).write();
    delete logger;
    fclose(fp);
    
    delete tree_writer;
    delete tree;
    delete staged_domain;
    delete domain_0;
    delete gen_data_0;

    printf("adan_kidney_vascularize ended.\n");
}

void adan_kidney_vascularize_cco_volumetric(double q0, vector<long long int> n_term,vector<AbstractConstraintFunction<double, int>*> gam,
  vector<AbstractConstraintFunction<double, int>*> eps_lim, vector<AbstractConstraintFunction<double, int>*> nu,
  vector<int> n_draw, vector<int> seed, vector<int> n_level_test, vector<int> n_terminal_trial, vector<double> d_lim_red_factor,
  vector<double> perfusion_area_factor, vector<double> close_neighborhood_factor, vector<double> mid_point_d_lim_factor, vector<int> n_bifurcation_test,
  vector<double> min_bif_angle, string cco_input, string suffix_id)
{
    printf("Beginning adan_kidney_vascularize.\n");

    // Timestamp   
    time_t now_time_t = time(nullptr);
    struct tm *now_tm = localtime(&now_time_t);
    char time_c_string[21];
    strftime(time_c_string, 20, "%d_%m_%Y_%H_%M_%S", now_tm);
    string time_string(time_c_string);

    // Folders location
    string input_folder {"/home/lfmc/HeMoLab/VascularizacaoRenal/Input/"};
    string output_folder {"/home/lfmc/HeMoLab/VascularizacaoRenal/Output/"};
    string steps_folder {"/home/lfmc/HeMoLab/VascularizacaoRenal/Steps/"};
    // Relevant input files
    //Geometry
    string input_vtk_pelvis {input_folder + "VTK/left_pelvis_06_07_2020_LFC_decimated.vtk"};
    string input_vtk_pyramids {input_folder + "VTK/left_pyramids_06_07_2020_LFC_decimated.vtk"};
    string input_vtk_core {input_folder + "VTK/left_core_06_07_2020_LFC_decimated.vtk"};
    string input_vtk_cortex {input_folder + "VTK/left_cortex_06_07_2020_LFC_decimated.vtk"};
    string input_vtk_surface {input_folder + "VTK/left_surface_06_07_2020_LFC_decimated.vtk"};
    //Vascular tree
    string input_cco {input_folder + "CCO/" + cco_input};
    // Output
    string output_base {"left_kidney_13_pyramids_" + suffix_id + "_" + "volumetric" + "_" + time_string};    
	
    AbstractVascularElement::VESSEL_FUNCTION function {AbstractVascularElement::VESSEL_FUNCTION::PERFORATOR};
    GeneratorData *gen_data_0 {new GeneratorData(n_level_test[0], n_terminal_trial[0], d_lim_red_factor[0], perfusion_area_factor[0], close_neighborhood_factor[0], mid_point_d_lim_factor[0], n_bifurcation_test[0], function, false)};
    DomainNVR *domain_0 {new DomainNVR(input_vtk_core, {input_vtk_pelvis, input_vtk_pyramids}, n_draw[0], seed[0], gen_data_0)};
    (*domain_0).setMinBifurcationAngle(min_bif_angle[0]);
    (*domain_0).setIsConvexDomain(false);

    StagedDomain *staged_domain {new StagedDomain()};
    (*staged_domain).addStage(n_term[0], domain_0);

    long long int n_term_total {n_term[0]};

    SingleVesselCCOOTree *tree {new SingleVesselCCOOTree(input_cco, gen_data_0, q0, gam[0], eps_lim[0], nu[0], 0.0, 1.0e-6)};
    tree->setIsInCm(true);
    StagedFRROTreeGenerator *tree_generator {new StagedFRROTreeGenerator(staged_domain, tree, n_term_total, gam, eps_lim, nu)};

    tree = {(SingleVesselCCOOTree *) (*tree_generator).resume(500, steps_folder)};
    make_terminal(tree);
    if(!tree) {
        printf("Failed to generate tree!\n");
        return;
    }
    (*tree).save(output_folder + "CCO/" + output_base + ".cco");
    VTKObjectTreeNodalWriter *tree_writer {new VTKObjectTreeNodalWriter()};
    (*tree_writer).write(output_folder + "VTP/" + output_base + ".vtp", tree);
    
    FILE* fp {fopen((output_folder + "LOG/" + output_base + ".log").c_str(), "w")};
    if(!fp) {
        fprintf(stderr, "Failed to create configuration log file.\n");
        exit(EXIT_FAILURE);
    }
    StagedFRROTreeGeneratorLogger *logger {new StagedFRROTreeGeneratorLogger(fp, tree_generator)};
    (*logger).write();
    delete logger;
    fclose(fp);

    delete tree_writer;
    delete tree;
    delete staged_domain;
    delete domain_0;
    delete gen_data_0;

    printf("adan_kidney_vascularize ended.\n");
}

int main(int argc, char *argv[])
{
    double q0 {14.};
    int n_level_test {16000};
    int n_terminal_trial {2000};
    double d_lim_red_factor {.95};
    double mid_point_d_lim_factor {.25};
    int n_bif_test {21};
    double vol_ref {100.0};
    double rad_ref {1.0};
    int factors[3][3] {{9990, 5, 5}, {9990, 10, 0}, {9990, 0, 10}};
    AbstractConstraintFunction<double,int> *gam {new ConstantConstraintFunction<double, int>(3.)};
    AbstractConstraintFunction<double,int> *eps_lim {new ConstantPiecewiseConstraintFunction<double, int>({0.4, 0.0}, {0, 7})};
    AbstractConstraintFunction<double,int> *nu {new ConstantConstraintFunction<double, int>(3.6)}; //cP

    vector<long long int> n_term_v({21});
    vector<int> n_draw_v(1, 2000);
    vector<int> n_level_test_v(1, n_level_test);
    vector<int> n_terminal_trial_v(1, n_terminal_trial);
    vector<double> d_lim_red_factor_v(1, d_lim_red_factor);
    vector<double> mid_point_d_lim_factor_v(1, mid_point_d_lim_factor);
    vector<int> n_bif_test_v(1, n_bif_test);
    vector<AbstractConstraintFunction<double,int> *> gam_v(1, gam);
    vector<AbstractConstraintFunction<double,int> *> eps_lim_v(1, eps_lim);
    vector<AbstractConstraintFunction<double,int> *> nu_v(1, nu);
    vector<double> min_bif_angle_v({(3. / 18.) * M_PI});
    vector<double> perfusion_area_factor_v(1, 1.0);
    vector<double> close_neighbourhood_factor_v({9.0});

    string input_cco{"renal_arteries_tri_simplified_04_02_2021.cco"};
    string suffix_id {"renal_artery_tri_radius"};

    vector<int> seed_v({1657});

    adan_kidney_vascularize_cco_functional(q0, n_term_v, gam_v, eps_lim_v, nu_v, n_draw_v, seed_v, n_level_test_v, n_terminal_trial_v, d_lim_red_factor_v,
					   perfusion_area_factor_v, close_neighbourhood_factor_v, mid_point_d_lim_factor_v, n_bif_test_v, min_bif_angle_v, vol_ref, rad_ref,
					   factors[0][0], factors[0][1], factors[0][2], input_cco, suffix_id);

    return 0;
}
