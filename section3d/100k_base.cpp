// STD Libs
#include<cmath>
#include<string>
#include<ctime>
#include<cstdio>
#include<cstdlib>
#include<vector>

// VItA Libs
#include<structures/domain/AbstractDomain.h>
#include<structures/domain/SimpleDomain.h>
#include<structures/domain/DomainNVR.h>
#include<structures/domain/PartiallyVascularizedDomain.h>
#include<structures/domain/StagedDomain.h>
#include<structures/domain/IntersectionVascularizedDomain.h>
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
#include<structures/domain/SimpleDomain2D.h>
#include<structures/vascularElements/AbstractVascularElement.h>

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

void vascularize_base(string output_base)
{
    string output_base_00 {output_base + "_p25"};
    string output_base_01 {output_base + "_p16"};
    string output_base_02 {output_base + "_p36"};

    int n_level_test {16000};
    int n_terminal_trial {2000};
    double d_lim_red_factor {.9};
    double mid_point_d_lim_factor {.25};
    double perfusion_area_factor {1.0};
    vector<double> close_neighborhood_factor({8.0, 4.0, 2.0});
    vector<int> seed({2049, 44, 57549});
    int n_bif_test {7};
    int n_draw {2000};
    vector<double> min_bif_angle(3, 3./18. * (M_PI));
    vector<long long int> n_term({833, 1666, 2501});
    // vector<long long int> n_term({10, 20, 70});
    AbstractConstraintFunction<double, int> *gam {new ConstantConstraintFunction<double, int>(3.0)};
    vector<AbstractConstraintFunction<double, int> *> gams(3, gam);
    AbstractConstraintFunction<double, int> *eps_lim {new ConstantPiecewiseConstraintFunction<double, int>({0.4, 0}, {0, 5})};
    vector<AbstractConstraintFunction<double, int> *> eps_lims(3, eps_lim);
    AbstractConstraintFunction<double,int> *nu {new ConstantConstraintFunction<double, int>(3.6)}; //cP
    vector<AbstractConstraintFunction<double, int> *> nus(3, nu);

    // Base geometry
    string input_vtk {input_vtk_folder + "pial_11_03_2021_hull.vtk"};

    AbstractVascularElement::VESSEL_FUNCTION function {AbstractVascularElement::VESSEL_FUNCTION::DISTRIBUTION};

    GeneratorData *gen_data_0 {new GeneratorData(n_level_test, n_terminal_trial, d_lim_red_factor, perfusion_area_factor, close_neighborhood_factor[0], mid_point_d_lim_factor, n_bif_test, function, true)};
    SimpleDomain *domain_0 {new SimpleDomain(input_vtk, n_draw, seed[0], gen_data_0)};
    (*domain_0).setMinBifurcationAngle(min_bif_angle[0]);
    (*domain_0).setIsConvexDomain(true);
    
    GeneratorData *gen_data_1 {new GeneratorData(n_level_test, n_terminal_trial, d_lim_red_factor, perfusion_area_factor, close_neighborhood_factor[1], mid_point_d_lim_factor, n_bif_test, function, false)};
    SimpleDomain *domain_1 {new SimpleDomain(input_vtk, n_draw, seed[1], gen_data_1)};
    (*domain_1).setMinBifurcationAngle(min_bif_angle[1]);
    (*domain_1).setIsConvexDomain(true);

    GeneratorData *gen_data_2 {new GeneratorData(n_level_test, n_terminal_trial, d_lim_red_factor, perfusion_area_factor, close_neighborhood_factor[2], mid_point_d_lim_factor, n_bif_test, function, false)};
    SimpleDomain *domain_2 {new SimpleDomain(input_vtk, n_draw, seed[2], gen_data_2)};
    (*domain_2).setMinBifurcationAngle(min_bif_angle[2]);
    (*domain_2).setIsConvexDomain(true);

    StagedDomain *staged_domain {new StagedDomain()};
    (*staged_domain).addStage(n_term[0], domain_0);
    (*staged_domain).addStage(n_term[1], domain_1);
    (*staged_domain).addStage(n_term[2], domain_2); 

    GeneratorData *gen_data_01 {new GeneratorData(n_level_test, n_terminal_trial, d_lim_red_factor, perfusion_area_factor, close_neighborhood_factor[2], mid_point_d_lim_factor, n_bif_test, function, false)};
    SimpleDomain *domain_01 {new SimpleDomain(input_vtk, n_draw, seed[2], gen_data_01)};
    (*domain_01).setMinBifurcationAngle(min_bif_angle[2]);
    (*domain_01).setIsConvexDomain(true);

    long long int total_term {0};
    for (auto it = n_term.begin(); it != n_term.end(); ++it) {
        total_term += *it;
    }

    long long int p16_surplus {8}, total_term_p16 {total_term + p16_surplus};
    StagedDomain *staged_domain_01 {new StagedDomain()};
    (*staged_domain_01).addStage(n_term[2] + p16_surplus, domain_01);
    (*staged_domain_01).setInitialStage(2);

    GeneratorData *gen_data_02 {new GeneratorData(n_level_test, n_terminal_trial, d_lim_red_factor, perfusion_area_factor, close_neighborhood_factor[2], mid_point_d_lim_factor, n_bif_test, function, false)};
    SimpleDomain *domain_02 {new SimpleDomain(input_vtk, n_draw, seed[2], gen_data_02)};
    (*domain_02).setMinBifurcationAngle(min_bif_angle[2]);
    (*domain_02).setIsConvexDomain(true);

    long long int p36_surplus {32 - p16_surplus}, total_term_p36 {total_term + p16_surplus + p36_surplus};
    StagedDomain *staged_domain_02 {new StagedDomain()};
    (*staged_domain_02).addStage(n_term[2] + p36_surplus, domain_02);
    (*staged_domain_02).setInitialStage(2);

    double q0 {2}; //mm3/s
    // diameter 150 micrometers 
    double r0 {0.075};
    point x0 {-2.49, -2.49, 1.15};
    
    StagedFRROTreeGenerator *tree_generator {new StagedFRROTreeGenerator(staged_domain, x0, r0, q0, total_term, gams, eps_lims, nus, 0., 1.e-5)};
    SingleVesselCCOOTree *tree = static_cast<SingleVesselCCOOTree *>((*tree_generator).getTree());
    (*tree).setIsInCm(false);
    (*tree_generator).generate(1, steps_base_folder);
    double dlim_last {(*tree_generator).getDLimLast()};

    (*tree).save(output_cco_folder + output_base_00 + ".cco");
    VTKObjectTreeNodalWriter *tree_writer {new VTKObjectTreeNodalWriter()};
    (*tree_writer).write(output_vtp_folder + output_base_00 + ".vtp", tree);    
    FILE* fp {fopen((output_log_folder + output_base_00 + ".log").c_str(), "w")};
    if(!fp) {
        fprintf(stderr, "Failed to create configuration log file.\n");
        exit(EXIT_FAILURE);
    }
    StagedFRROTreeGeneratorLogger *logger {new StagedFRROTreeGeneratorLogger(fp, tree_generator)};
    (*logger).write();
    delete logger;
    fclose(fp);


    StagedFRROTreeGenerator *tree_generator_01 {new StagedFRROTreeGenerator(staged_domain_01, tree, total_term_p16, gams, eps_lims, nus)};
    (*tree_generator_01).setDLim(dlim_last);
    SingleVesselCCOOTree *tree_01 = static_cast<SingleVesselCCOOTree *>((*tree_generator_01).getTree());
    (*tree_01).setIsInCm(false);
    (*tree_generator_01).resume(1, steps_base_folder);
    dlim_last = (*tree_generator_01).getDLimLast();

    (*tree_01).save(output_cco_folder + output_base_01 + ".cco");
    (*tree_writer).write(output_vtp_folder + output_base_01 + ".vtp", tree_01);
    FILE* fp_01 {fopen((output_log_folder + output_base_01 + ".log").c_str(), "w")};
    if(!fp_01) {
        fprintf(stderr, "Failed to create configuration log file.\n");
        exit(EXIT_FAILURE);
    }
    StagedFRROTreeGeneratorLogger *logger_01 {new StagedFRROTreeGeneratorLogger(fp_01, tree_generator_01)};
    (*logger_01).write();
    delete logger_01;
    fclose(fp_01);

    StagedFRROTreeGenerator *tree_generator_02 {new StagedFRROTreeGenerator(staged_domain_02, tree_01, total_term_p36, gams, eps_lims, nus)};
    (*tree_generator_02).setDLim(dlim_last);
    SingleVesselCCOOTree *tree_02 = static_cast<SingleVesselCCOOTree *>((*tree_generator_02).getTree());
    (*tree_02).setIsInCm(false);
    (*tree_generator_02).resume(1, steps_base_folder);
    
    (*tree).save(output_cco_folder + output_base_02 + ".cco");
    (*tree_writer).write(output_vtp_folder + output_base_02 + ".vtp", tree_02);
    FILE* fp_02 {fopen((output_log_folder + output_base_02 + ".log").c_str(), "w")};
    if(!fp_02) {
        fprintf(stderr, "Failed to create configuration log file.\n");
        exit(EXIT_FAILURE);
    }
    StagedFRROTreeGeneratorLogger *logger_02 {new StagedFRROTreeGeneratorLogger(fp, tree_generator_02)};
    (*logger_02).write();
    delete logger_02;
    fclose(fp_02);
      
    delete tree_writer;
    delete tree_generator;
    delete tree_generator_01;
    delete tree_generator_02;
    delete staged_domain_02;
    delete domain_02;
    delete gen_data_02;
    delete staged_domain_01;
    delete domain_01;
    delete gen_data_01;
    delete staged_domain;
    delete domain_2;
    delete gen_data_2;
    delete domain_1;
    delete gen_data_1;
    delete domain_0;
    delete gen_data_0;
    delete gam;
    delete nu;
    delete eps_lim;
}

int main(int argc, char *argv[])
{
    vascularize_base("100k_base");
    return 0;
} 
