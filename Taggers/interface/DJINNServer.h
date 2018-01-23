
#include <vector>
#include "DNN/TensorFlow/interface/TensorFlow.h"

using namespace std;

class DJINNServer{

    public:
        DJINNServer();
        DJINNServer(const string model_dir,
                    unsigned pixels,
                    unsigned channels,
                    unsigned eng_features);

        ~DJINNServer();

        float evaluate(const vector<vector<vector<float>>> &image,
                        const vector<float> &eng_features);

    private:

        unsigned pixels_;
        unsigned channels_;
        unsigned eng_features_;

/*
        dnn::tf::Graph* g_;
        dnn::tf::Tensor* x_im_;
        dnn::tf::Tensor* x_ef_;

        dnn::tf::Tensor* kp_conv_;
        dnn::tf::Tensor* kp_hidd_;
        dnn::tf::Tensor* inference_;

        dnn::tf::Tensor* y_;
*/
        tensorflow::MetaGraphDef* metaGraph_;
        tensorflow::Session* session_;
};


DJINNServer::DJINNServer(const string model_dir, unsigned pixels,
                         unsigned channels, unsigned eng_features){
    
    //Input properties
    pixels_  = pixels;
    channels_ = channels;
    eng_features_ = eng_features;

    metaGraph_ = tensorflow::loadMetaGraph(model_dir);
    session_ = tensorflow::createSession(metaGraph_,model_dir);
/*    
    //Tensorflow model setup
    //Load computational graph from model file
    g_ = new dnn::tf::Graph(model_dir);

    //Define input tensors and their shapes
    dnn::tf::Shape x_im_Shape[] = {1,pixels_,pixels_,channels_};
    dnn::tf::Shape x_ef_Shape[] = {1,eng_features_};
    dnn::tf::Shape kp_Shape[]   = {1};

    x_im_ = g_->defineInput(new dnn::tf::Tensor("im_in:0", 4, x_im_Shape));
    x_ef_ = g_->defineInput(new dnn::tf::Tensor("ef_in:0", 2, x_ef_Shape));

    kp_conv_ = g_->defineInput(new dnn::tf::Tensor("kp_conv:0", 1, kp_Shape));
    kp_hidd_ = g_->defineInput(new dnn::tf::Tensor("kp_hidd:0", 1, kp_Shape));
    inference_ = g_->defineInput(new dnn::tf::Tensor("inference:0", 1, kp_Shape));

    //Define output tensor
    y_ = g_->defineOutput(new dnn::tf::Tensor("y_prob:0"));
*/


}

float DJINNServer::evaluate(const vector<vector<vector<float>>> &image,
                        const vector<float> &eng_features){
/*
    //Settings for inference
    kp_conv_->setValue<float>(0,1.0);
    kp_hidd_->setValue<float>(0,1.0);
    inference_->setValue<bool>(0,true);

    //Fill input image placeholder
    for (unsigned i=0;i<pixels_;i++){
        for (unsigned j=0;j<pixels_;j++){
            for (unsigned k=0;k<channels_;k++){
                x_im_->setValue<float>(0,i,j,k,image[i][j][k]);
            }
        }
    }

    //Fill engineered features placeholder
    for (unsigned i=0;i<eng_features_;i++){
        x_ef_->setValue<float>(0,i,eng_features[i]);
    }

    //Execute computational graph
    g_->eval();

    //Retrieve value and return
    return y_->getValue<float>(0,1);
*/

    return 1.0;

}



