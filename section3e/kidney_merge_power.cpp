#include<cstdio>
#include<cstdlib>
#include<cstring>
#include<string>
#include<vector>

#include<core/TreeMerger.h>
#include<structures/tree/SingleVesselCCOOTree.h>
#include<core/GeneratorData.h>
#include<constrains/AbstractConstraintFunction.h>
#include<constrains/ConstantConstraintFunction.h>
#include<constrains/ConstantPiecewiseConstraintFunction.h>
#include<io/VTKObjectTreeNodalWriter.h>
#include<structures/tree/VolumetricCostEstimator.h>
#include<structures/tree/DoublePowerCostEstimator.h>
#include<structures/tree/AdimSproutingVolumetricCostEstimator.h>
#include<structures/tree/SproutingVolumetricCostEstimator.h>
#include<structures/tree/PowerCostEstimator.h>
#include<structures/tree/LinearCombinationCostEstimator.h>

void merge(string input_cco, vector<string> parts, AbstractConstraintFunction<double, int> *gam, string output_base) 
{
    AbstractConstraintFunction<double, int> *eps_lim {new ConstantPiecewiseConstraintFunction<double, int>({0.4, 0}, {0, 7})};
    AbstractConstraintFunction<double, int> *nu {new ConstantConstraintFunction<double, int>(3.6)};
    GeneratorData *gen_data {new GeneratorData()};

    SingleVesselCCOOTree *tree {new SingleVesselCCOOTree(input_cco, gen_data, gam, eps_lim, nu)};
    (*tree).setIsFL(true);
    (*tree).setIsInCm(true);

    TreeMerger *merger {new TreeMerger(tree, parts)};
    (*merger).mergeFast();
    VTKObjectTreeNodalWriter *writer {new VTKObjectTreeNodalWriter()};
    (*tree).save(output_base + ".cco");
    (*writer).write(output_base + ".vtp", tree);
    
    delete writer;
    delete merger;
    delete tree;
    delete gen_data;
    delete nu;
    delete eps_lim;
}

int main(int argc, char *argv[])
{   
    char cco_path[500];
    char points_path[500];
    char output_path[500];
    char suffix_c[100];
    double gamma_value;
    
    AbstractConstraintFunction<double, int> *gam {nullptr};

    FILE *fp = fopen(argv[1], "r");
    if(!fp) {
        fprintf(stderr, "Failed to open seed file.\n");
        exit(EXIT_FAILURE);
    }
    fscanf(fp, "%s\n", cco_path);
    fscanf(fp, "%s\n", points_path);
    fscanf(fp, "%s\n", output_path);
    fscanf(fp, "%lf\n", &gamma_value);
    gam = new ConstantConstraintFunction<double, int>(gamma_value);
    fscanf(fp, "%s\n", suffix_c);
    printf("%s\n", suffix_c);
    fclose(fp);

    string input_cco {string(cco_path) + "kidney_base_" + string(suffix_c) + ".cco"};
    string output_base {string(output_path) + "kidney_merged_" + string(suffix_c)};
    vector<string> input_points;

    for (int i = 1; i<=8; ++i) {
        char no[5];
        sprintf(no, "%02d", i);
        input_points.push_back(string(points_path) + "kidney_region_" + string(no) + "_" + string(suffix_c) + ".points"); 
    }
    
    merge(input_cco, input_points, gam, output_base);

    delete gam;

    return 0;
}