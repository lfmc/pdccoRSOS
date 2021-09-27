#include<cstdio>
#include<cstring>
#include<string>

#include<structures/vascularElements/AbstractVascularElement.h>
#include<structures/tree/SingleVesselCCOOTree.h>
#include<constrains/ConstantConstraintFunction.h>
#include<constrains/ConstantPiecewiseConstraintFunction.h>
#include<structures/domain/SimpleDomain.h>
#include<structures/domain/StagedDomain.h>
#include<core/GeneratorData.h>
#include<core/TreeMerger.h>
#include<core/StagedFRROTreeGenerator.h>
#include<io/VTKObjectTreeNodalWriter.h>
#include<io/StagedFRROTreeGeneratorLogger.h>

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

string vascularize_part(string output_base, string input_cco, string input_vtk, long long int base_terms, long long int region_terms, double dlim)
{
    int n_stages {7};
    int n_level_test {16000};
    int n_terminal_trial {2000};
    double d_lim_red_factor {.9};
    double mid_point_d_lim_factor {.25};
    double perfusion_area_factor {1.0};
    vector<double> close_neighborhood_factor({1.0, 0.5, 0.25, 0.125});
    vector<int> seed({2049, 44, 57549, 1917});
    int n_bifurcation_test {7};
    int n_draw {2000};
    double min_bif_angle {(1./6.) * (M_PI)};
    AbstractConstraintFunction<double, int> *gam {new ConstantConstraintFunction<double, int>(3.0)};
    vector<AbstractConstraintFunction<double, int> *> gams(n_stages, gam);
    AbstractConstraintFunction<double, int> *eps_lim {new ConstantPiecewiseConstraintFunction<double, int>({0.4, 0}, {0, 5})};
    vector<AbstractConstraintFunction<double, int> *> eps_lims(n_stages, eps_lim);
    AbstractConstraintFunction<double,int> *nu {new ConstantConstraintFunction<double, int>(3.6)}; //cP
    vector<AbstractConstraintFunction<double, int> *> nus(n_stages, nu);
    AbstractVascularElement::VESSEL_FUNCTION function {AbstractVascularElement::VESSEL_FUNCTION::DISTRIBUTION};

    // long long base_terms {5000};
    lldiv_t term_division = lldiv(region_terms, 10);
    vector<long long int> n_term({term_division.quot, 2 * term_division.quot, 3 * term_division.quot,  (4 * term_division.quot) + term_division.rem});

    GeneratorData *gen_data_0 {new GeneratorData(n_level_test, n_terminal_trial, d_lim_red_factor, perfusion_area_factor, close_neighborhood_factor[0], mid_point_d_lim_factor, n_bifurcation_test, function, false)};
    SimpleDomain *domain_0 {new SimpleDomain(input_vtk, n_draw, seed[0], gen_data_0)};
    (*domain_0).setMinBifurcationAngle(min_bif_angle);
    (*domain_0).setIsConvexDomain(false);

    GeneratorData *gen_data_1 {new GeneratorData(n_level_test, n_terminal_trial, d_lim_red_factor, perfusion_area_factor, close_neighborhood_factor[1], mid_point_d_lim_factor, n_bifurcation_test, function, false)};
    SimpleDomain *domain_1 {new SimpleDomain(input_vtk, n_draw, seed[1], gen_data_1)};
    (*domain_1).setMinBifurcationAngle(min_bif_angle);
    (*domain_1).setIsConvexDomain(false);

    GeneratorData *gen_data_2 {new GeneratorData(n_level_test, n_terminal_trial, d_lim_red_factor, perfusion_area_factor, close_neighborhood_factor[2], mid_point_d_lim_factor, n_bifurcation_test, function, false)};
    SimpleDomain *domain_2 {new SimpleDomain(input_vtk, n_draw, seed[2], gen_data_2)};
    (*domain_2).setMinBifurcationAngle(min_bif_angle);
    (*domain_2).setIsConvexDomain(false);

    GeneratorData *gen_data_3 {new GeneratorData(n_level_test, n_terminal_trial, d_lim_red_factor, perfusion_area_factor, close_neighborhood_factor[3], mid_point_d_lim_factor, n_bifurcation_test, function, false)};
    SimpleDomain *domain_3 {new SimpleDomain(input_vtk, n_draw, seed[3], gen_data_3)};
    (*domain_3).setMinBifurcationAngle(min_bif_angle);
    (*domain_3).setIsConvexDomain(false);

    StagedDomain *staged_domain {new StagedDomain()};
    (*staged_domain).addStage(n_term[0], domain_0);
    (*staged_domain).addStage(n_term[1], domain_1);
    (*staged_domain).addStage(n_term[2], domain_2);
    (*staged_domain).addStage(n_term[3], domain_3);
    int mergeStage {3};
    (*staged_domain).setInitialStage(mergeStage);

    printf("Trying to read cco file.\n");
    SingleVesselCCOOTree *tree {new SingleVesselCCOOTree(input_cco, gen_data_0, gam, eps_lim, nu)};
    (*tree).setIsInCm(false);
    StagedFRROTreeGenerator *tree_generator {new StagedFRROTreeGenerator(staged_domain, tree, base_terms + region_terms, gams, eps_lims, nus)};
    (*tree_generator).setDLim(dlim);
       
    string output_points = output_points_folder + output_base + ".points";
    FILE *fp_points = fopen(output_points.c_str(), "wb");
    if (!fp_points) {
        fprintf(stderr, "Failed to create points files!\n");
        exit(EXIT_FAILURE);
    }
    tree = {(SingleVesselCCOOTree *) (*tree_generator).resumeSavePointsMidPoint(500, steps_base_folder, fp_points)};
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
    delete gam;
    delete eps_lim;
    delete nu;

    return output_points;
}


int main(int argc, char *argv[]) 
{
    if (argc != 3) {
        fprintf(stderr, "Not enough arguments.\n");
        exit(EXIT_FAILURE);
    }

    string base_folder {"/scratch/hmlffr2/luis.cury/"};
    double dlim {0};
    long long int n_part_term {0};
    string input_cco {input_cco_folder + "100k_base_" + string(argv[1]) + ".cco"};
    string input_vtk {input_vtk_folder + "network_"};
    int part {0};
    sscanf(argv[2], "%d", &part);
    char which_part[5];
    sprintf(which_part, "%03d", part);
    input_vtk += string(argv[1]) + "_" + string(which_part) + ".vtk";
    string output_base {"100k_" + string(argv[1]) + "_" + string(which_part)};
    long long int base_terms {5000};

    if (!strcmp(argv[1], "p16")) {
        dlim = 0.146133;
        n_part_term = 5937;
        base_terms += 8;
    }
    else if (!strcmp(argv[1], "p25")) {
        dlim = 0.146211;
        n_part_term = 3800;
    }
    else if (!strcmp(argv[1], "p36")) {
        dlim = 0.145900;
        n_part_term = 2638;
        base_terms += 32;
    }
    else {
        fprintf(stderr, "Invalid case.\n");
    }
    printf("output_base = %s\n", output_base.c_str());
    printf("input_cco = %s\n", input_cco.c_str());
    printf("input_vtk = %s\n", input_vtk.c_str());
    vascularize_part(output_base, input_cco, input_vtk, base_terms, n_part_term, dlim);
    return 0;
}
