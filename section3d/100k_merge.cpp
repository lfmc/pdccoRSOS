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

void merge(string input_cco, vector<string> parts, string output_base) 
{
    
    
    GeneratorData *gen_data {new GeneratorData()};
    AbstractConstraintFunction<double, int> *gam {new ConstantConstraintFunction<double, int>(3.0)};
    AbstractConstraintFunction<double, int> *eps_lim {new ConstantPiecewiseConstraintFunction<double, int>({0.4, 0}, {0, 5})};
    AbstractConstraintFunction<double,int> *nu {new ConstantConstraintFunction<double, int>(3.6)}; //cP

    SingleVesselCCOOTree *tree {new SingleVesselCCOOTree(input_cco, gen_data, gam, eps_lim, nu)};
    (*tree).setIsFL(false);
    (*tree).setIsInCm(false);
    
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
    if (argc != 3) {
        fprintf(stderr, "Usage: ./exec location_folder part(int)\n");
        exit(EXIT_FAILURE);
    }
    int type {0};
    string folder_base {argv[1]};
    sscanf(argv[2], "%d", &type);
    printf("p%02d\n", type);
    char type_c[10];
    sprintf(type_c, "p%02d", type);
    string input_cco {folder_base + "/100k_base_" + string(type_c) + ".cco"};
    vector<string> input_points;
    for (int i = 0; i < type; ++i){
        char file_suffix[30];
        sprintf(file_suffix, "p%02d_%03d", type, i+1);
        string points_input {folder_base + "/100k_" + string(file_suffix) + ".points"};
        printf("%s\n", points_input.c_str());
        string points_file {points_input};
        input_points.push_back(points_file);
    }
    string output_base {folder_base + "/100k_merged_" + string(type_c)};
    merge(input_cco, input_points, output_base);
    return 0;
}
