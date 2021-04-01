#include <iostream>
#include <vector>
#include <string>

#include <TFile.h>
#include <TTree.h>

#include "SolidBox.hh"
#include "GeometrySimple.hh"
#include "ComponentConstant.hh"
#include "MediumMagboltz.hh"
#include "Sensor.hh"
#include "TrackSrim.hh"
#include "Random.hh"
#include "Plotting.hh"
#include "AvalancheMC.hh"
#include "TrackSrim.hh"
#include "Random.hh"

#include "TrackTrimSQLite.hpp"

int main()
{

    const double world_size = 100; //(cm)

    auto gas = new Garfield::MediumMagboltz();
    auto box = new Garfield::SolidBox(0, 0, 0, world_size, world_size, world_size);
    auto geo = new Garfield::GeometrySimple();
    auto comp = new Garfield::ComponentConstant();
    auto sensor = new Garfield::Sensor();

    geo->AddSolid(box, gas);
    comp->SetGeometry(geo);
    sensor->AddComponent(comp);

    //Gasfile
    const std::string gasfile = "../input/He(90)+CH4(10)_1000.gas";
    gas->LoadGasFile(gasfile);

    // Gas property
    double a_eff, z_eff, w_eff;
    a_eff = 4 * 0.9 + 16 * 0.1;
    z_eff = 2 * 0.9 + 10 * 0.1;
    w_eff = (2 * 41.3 * 0.9 + 10 * 30 * 0.1) / z_eff;

    // TrackClass
    auto track = new GarfieldSuppl::TrackTrimSQLite();
    track->ReadFile("./hoge.sqlite");
    track->SetSensor(sensor);
    track->SetWorkFunction(w_eff);
    track->SetTargetClusterSize(1);
    track->SetTransferProbability(0.1);
    // track->EnableDebugging();

    //value for Tree
    std::vector<double> vx, vy, vz, vt;
    std::vector<int> vn;
    std::vector<double> ve_cls, ve_ion;

    TFile fOut("tracktrim_new.root", "recreate");
    TTree tr("tr", "");
    tr.Branch("x", &vx);
    tr.Branch("y", &vy);
    tr.Branch("z", &vz);
    tr.Branch("t", &vt);
    tr.Branch("n", &vn);
    tr.Branch("e_cls", &ve_cls);
    tr.Branch("e_ion", &ve_ion);

    for (int iTrack = 0; iTrack < 1000; ++iTrack)
    {
        track->SetKineticEnergy(1.0e+6);
        std::cout << " " << iTrack << std::endl;
        track->NewTrack(0, 0, 0, 0, 1, 0, 0);

        double xc, yc, zc, tc;
        int ne_c;
        double e_cls, e_ion;

        vx.clear();
        vy.clear();
        vz.clear();
        vt.clear();
        vn.clear();
        ve_cls.clear();
        ve_ion.clear();

        while (track->GetCluster(xc, yc, zc, tc, ne_c, e_cls, e_ion))
        {
            vx.push_back(xc);
            vy.push_back(yc);
            vz.push_back(zc);
            vt.push_back(tc);
            vn.push_back(ne_c);
            ve_cls.push_back(e_cls);
            ve_ion.push_back(e_ion);
        }

        tr.Fill();
    }

    tr.Write();

    // delete all of Garfield objects
    delete gas;
    delete box;
    delete geo;
    delete comp;
    delete sensor;
    delete track;

    return 0;
}