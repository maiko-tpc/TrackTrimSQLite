#pragma once

#include <memory>
#include "Track.hh"

#include "TrackGenerator.hpp"

namespace GarfieldSuppl
{

    using namespace Garfield;

    class TrackTrimSQLite : public Track
    {

    public:
        /// Constructor
        TrackTrimSQLite();
        /// Destructor
        virtual ~TrackTrimSQLite();

        /// Set/get the W value [eV].
        void SetWorkFunction(const double w) { m_work = w; }
        double GetWorkFunction() const { return m_work; }

        double GetTransferProbability() const
        {
            return m_generator->GetTransferProbability();
        };

        void SetTransferProbability(const double prob)
        {
            m_generator->SetTransferProbability(prob);
        };

        double GetEnergyMarginRatio() const
        {
            return m_generator->GetEnergyMarginRatio();
        };
        void SetEnergyMarginRatio(const double ratio_to_de)
        {
            return m_generator->SetEnergyMarginRatio(ratio_to_de);
        };

        void SetTargetClusterSize(const int n) { m_nsize = n; }
        int GetTargetClusterSize() const { return m_nsize; }

        void SetClustersMaximum(const int n) { m_maxclusters = n; }
        int GetClustersMaximum() const { return m_maxclusters; }

        bool ReadFile(const std::string &file);

        virtual bool NewTrack(const double x0, const double y0, const double z0,
                              const double t0, const double dx0, const double dy0,
                              const double dz0);
        virtual bool GetCluster(double &xcls, double &ycls, double &zcls,
                                double &tcls, int &n, double &e, double &extra);

    protected:
        /// Work function [eV]
        double m_work = -1.;

        /// Maximum number of clusters allowed (infinite if 0)
        int m_maxclusters = -1;

        /// Index of the next cluster to be returned
        unsigned int m_currcluster;

        /// Targeted cluster size
        int m_nsize = -1;

        struct cluster
        {
            double x, y, z, t; // Cluster location and time
            double ec;         // Energy spent to make the clusterec
            double kinetic;    // Ion energy when cluster was created
            int electrons;     // Number of electrons in this cluster
        };
        std::vector<cluster> m_clusters;

        std::unique_ptr<TrackGeneratorC> m_generator;

        bool NewClusterPushable() const
        {
            return m_maxclusters < 0 || m_clusters.size() < m_maxclusters;
        }

        bool IsInside(const double x, const double y, const double z);
    };
}