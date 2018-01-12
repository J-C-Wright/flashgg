

#include <iostream>
#include <vector>
#include <math.h>

void print_constituent_vector(std::vector<std::vector<float>> &reshaped){

    for (unsigned j=0;j<reshaped.size();j++){
        for (unsigned k=0;k<4;k++){
            std::cout << reshaped[j][k] << " ";
        }
        std::cout << std::endl;
    }
};

void constituent_vector_from_flat(std::vector<float> &flat_vec,
                                  std::vector<std::vector<float>> &reshaped){


    unsigned n_cons = flat_vec.size()/4;
    for (unsigned j=0;j<n_cons;j++){

        std::vector<float>::const_iterator first = flat_vec.begin() + j*4;
        std::vector<float>::const_iterator last  = flat_vec.begin() + (j+1)*4;
        std::vector<float> subvec(first,last);

        reshaped.push_back(subvec);
    }

};
void constituent_vector_from_tree(unsigned entry, TTree* tree, TString branch_name, 
                                  std::vector<std::vector<float>> &reshaped){

    std::vector<float> *cons = new std::vector<float>();
    tree->SetBranchAddress(branch_name,&cons); 

    tree->GetEntry(entry);

    unsigned n_cons = cons->size()/4;

    for (unsigned j=0;j<n_cons;j++){

        std::vector<float>::const_iterator first = cons->begin() + j*4;
        std::vector<float>::const_iterator last  = cons->begin() + (j+1)*4;
        std::vector<float> subvec(first,last);

        reshaped.push_back(subvec);
    }

};

void image_from_vector(std::vector<std::vector<std::vector<float>>> &image,
                       std::vector<std::vector<float>> &reshaped,
                       unsigned n_pixels,float radius, unsigned merge){

    float cos_d = cos(M_PI/float(n_pixels));
    float sin_d = sin(M_PI/float(n_pixels));

    float n_const  = 0.0;
    float total_pt = 0.0;


    for (unsigned i=0;i<reshaped.size();i++){

        float R = sqrt(pow(reshaped[i][0],2)+pow(reshaped[i][1],2));  
        float phi = atan2(sin_d*reshaped[i][0]+cos_d*reshaped[i][1],
                           cos_d*reshaped[i][0]-sin_d*reshaped[i][1]);

        unsigned i_R   = unsigned(n_pixels*R/radius);
        unsigned i_phi = unsigned(n_pixels*(phi+M_PI)/(2*M_PI));
        unsigned i_c = abs(abs(reshaped[i][2])-1);

        if (i_R >= n_pixels) continue;
        if (i_phi >= n_pixels) continue;

        image[i_phi][i_R][i_c] += reshaped[i][3];
        image[i_phi][i_R][2] += 1.0;

        total_pt += reshaped[i][3];
        n_const += 1.0;

    }

    //pt channel norm
    for (unsigned i=0;i<n_pixels;i++){
        for (unsigned j=0;j<n_pixels;j++){
            for (unsigned c=0;c<2;c++){
                image[i][j][c] /= total_pt;
            }
        }
    }

    //multiplicity norm
    for (unsigned i=0;i<n_pixels;i++){
        for (unsigned j=0;j<n_pixels;j++){
            image[i][j][2] /= n_const;
        }
    }

    //Central merging
    if (merge > 0) {
        unsigned central_bins = n_pixels/merge;
        for (unsigned c=0;c<3;c++){
            for (unsigned i=0;i<central_bins;i++){
                float sum = 0.0;
                for (unsigned j=0;j<merge;j++){
                    sum += image[i*merge+j][0][c];
                }
                for (unsigned j=0;j<merge;j++){
                    image[i*merge+j][0][c] = sum;
                }
            }
        }
    }

    
};



std::vector<std::vector<std::vector<float>>> blank_image(unsigned n_pixels, unsigned n_channels){

    std::vector<std::vector<std::vector<float>>> image;
    std::vector<float> pixel(n_channels,0.0);

    for (unsigned i=0;i<n_pixels;i++){

        std::vector<std::vector<float>> row;

        for (unsigned j=0;j<n_pixels;j++){
            row.push_back(pixel);
        }
        image.push_back(row);
    }

    return image;
};

void print_image(std::vector<std::vector<std::vector<float>>> image){

    unsigned n_pixels = image.size();
    unsigned n_chans = image[0][0].size();

    for (unsigned c=0;c<n_chans;c++){
        for (unsigned i=0;i<n_pixels;i++){
            for (unsigned j=0;j<n_pixels;j++){
                std::cout << image[i][j][c] << " ";
            }
            std::cout << std::endl;
        }
        std::cout << std::endl;
    }
};

/*
void test(){


    TString path = "/Users/Jack/Dropbox/Code/DJINN/dataset/VBFHToGG_M125_13TeV_amcatnlo_pythia8.root";
    TString selection = "( leadJetPt > 30 && subleadJetPt > 20 && lead_pt_m > (1/4.0) && sublead_pt_m > (1/5.0) && (leadJetEta < 4.7 && leadJetEta > -4.7) && (subleadJetEta < 4.7 && subleadJetEta > -4.7) && mass_jj > 100 && mass_gg > 100 && mass_gg < 180 )";
    TFile *file = TFile::Open(path);

    TString tree_path = "flashggDJINNTreeMaker/JetData";
    TTree *tree_full = (TTree*)file->Get(tree_path);

    TTree* tree = tree_full->CopyTree(selection);

    unsigned n_pixels = 8;
    unsigned n_channels = 3;
    TString branch_name = "leadConstituents";

    std::vector<std::vector<std::vector<float>>> mean_image = blank_image(n_pixels,n_channels);

    unsigned n_events = 1000;
    for (unsigned i=0;i<n_events;i++){

        std::vector<std::vector<float>> reshaped;
        constituent_vector_from_tree(i,tree,branch_name,reshaped);

        std::vector<std::vector<std::vector<float>>> image = blank_image(n_pixels,n_channels);
        image_from_vector(image,reshaped,n_pixels,0.4,2);

        for (unsigned i=0;i<n_pixels;i++){
            for (unsigned j=0;j<n_pixels;j++){
                for (unsigned c=0;c<3;c++){
                    mean_image[i][j][c] += image[i][j][c];
                }
            }
        }
    }

    for (unsigned i=0;i<mean_image.size();i++){
        for (unsigned j=0;j<mean_image[i].size();j++){
            for (unsigned k=0;k<mean_image[i][j].size();k++){
                mean_image[i][j][k] /= float(n_events);
            }
        }
    }

    print_image(mean_image);

}
*/


