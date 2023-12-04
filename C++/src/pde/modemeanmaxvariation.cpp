#include "modemeanmaxvariation.h"

ModeMeanMaxVariation::ModeMeanMaxVariation() {}

/*
    Array<float> mean;
    shared_ptr<BaseVector> mean_vec = CreateBaseVector(evecs[0]->Size(), true, 1);
    mean.SetSize(evecs[0]->Size());

    for(int i = 0; i < mean.Size(); i++)
    {
        mean[i] = 0.0f;
        for(int j = 0; j < evecs.Size(); j++)
        {
            mean[i] += std::abs(evecs[j]->FV<std::complex<double>>()[i]);
        }
        mean[i] /= evecs.Size();
    }

    std::vector<float> var;
    for(int i = 0; i < mean.Size(); i++)
    {
        float s = 0.0f;
        for(int j = 0; j < evecs.Size(); j++)
        {
             float tmp = std::abs(evecs[j]->FV<std::complex<double> >()[i]) - mean[i];
             tmp *= tmp;
             s += tmp;
        }
        s /= evecs.Size();
        var.push_back(std::sqrt(s));
        mean_vec->FV<std::complex<double>>()[i] = std::sqrt(s);
        std::cout << std::sqrt(s) << std::endl;
    }

    Array<shared_ptr<BaseVector>> arr_mean;
    arr_mean.Append(mean_vec);

    SolutionMapper sm(ma);
    //sm.mapSolutionToMesh(elasticityFes.getFes(), arr_mean);
    sm.mapSolutionToMesh(elasticityFes.getFes(), ess.getEvecs());

    //mean_vec->Print(std::cout);
    std::cout << " - > " << mean_vec->Size() << ", " << evecs[0]->Size() << std::endl;
    */
