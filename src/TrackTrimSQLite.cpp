#include <GarfieldConstants.hh>
#include <Sensor.hh>
#include <Random.hh>

#include "TrackTrimSQLite.hpp"

namespace GarfieldSuppl
{

    using namespace Garfield;

    /// Constructor
    TrackTrimSQLite::TrackTrimSQLite()
        : Track(), m_generator(new TrackGeneratorC())
    {
        m_className = "TrackTrimSQLite";
        m_generator->SetRandomGenerator(RndmUniform);
    }

    /// Destructor
    TrackTrimSQLite::~TrackTrimSQLite() {}

    bool TrackTrimSQLite::ReadFile(const std::string &file)
    {
        return m_generator->SetFileName(file);
    }

    bool TrackTrimSQLite::NewTrack(const double x0, const double y0, const double z0,
                                   const double t0, const double dx0, const double dy0,
                                   const double dz0)
    {

        // Reset the cluster count
        m_currcluster = 0;
        m_clusters.clear();

        const std::string hdr = m_className + "::NewTrack: ";

        if (!m_generator->IsAccesible())
        {
            std::cerr << hdr << "\n    SQLite collision DB is inaccessible.\n";
            return false;
        }

        // Verify that a sensor has been set.
        if (!m_sensor)
        {
            std::cerr << hdr << "\n    Sensor is not defined.\n";
            return false;
        }

        // Get the bounding box.
        double xmin = 0., ymin = 0., zmin = 0.;
        double xmax = 0., ymax = 0., zmax = 0.;
        if (!m_sensor->GetArea(xmin, ymin, zmin, xmax, ymax, zmax))
        {
            std::cerr << hdr << "\n    Drift area is not set.\n";
            return false;
        }
        else if (x0 < xmin || x0 > xmax ||
                 y0 < ymin || y0 > ymax || z0 < zmin || z0 > zmax)
        {
            std::cerr << hdr << "\n    Initial position outside bounding box.\n";
            return false;
        }

        // Make sure the initial position is inside an ionisable medium.
        Medium *medium = NULL;
        if (!m_sensor->GetMedium(x0, y0, z0, medium))
        {
            std::cerr << hdr << "\n    No medium at initial position.\n";
            return false;
        }
        else if (!medium->IsIonisable())
        {
            std::cerr << hdr << "\n    Medium at initial position is not ionisable.\n";
            return false;
        }

        // Normalise and store the direction.
        const double normdir = sqrt(dx0 * dx0 + dy0 * dy0 + dz0 * dz0);
        double xdir = dx0;
        double ydir = dy0;
        double zdir = dz0;
        if (normdir < Small)
        {
            if (m_debug)
            {
                std::cout << hdr << "\n    Direction vector has zero norm.\n"
                          << "    Initial direction is randomized.\n";
            }
            // Null vector. Sample the direction isotropically.
            RndmDirection(xdir, ydir, zdir);
        }
        else
        {
            // Normalise the direction vector.
            xdir /= normdir;
            ydir /= normdir;
            zdir /= normdir;
        }

        // Make sure all necessary parameters have been set.
        if (m_energy < Small)
        {
            std::cerr << hdr << "\n    Initial particle energy not set.\n";
            return false;
        }
        else if (m_work < Small)
        {
            std::cerr << hdr << "\n    Work function not set.\n";
            return false;
        }
        // Check the initial energy.
        const double ekin0 = GetKineticEnergy();
        if (ekin0 < m_work)
        {
            if (m_debug)
            {
                std::cout << hdr << "Initial kinetic energy E = " << ekin0
                          << " eV E < W; particle stopped.\n";
            }
            return true;
        }

        // // Get an upper limit for the track length.
        // const double tracklength = 10 * Interpolate(ekin0, m_ekin, m_range);

        // Header of debugging output.
        if (m_debug)
        {
            std::cout << hdr << "Track generation with the following parameters:\n";
            // const unsigned int nTable = m_ekin.size();
            printf("      DB file              %s\n", m_generator->GetFileName().c_str());
            printf("      Particle kin. energy %g keV\n", ekin0);
            //printf("      Particle mass        %g MeV\n", 1.e-6 * m_mass);
            printf("      Particle charge      %g\n", m_q);
            printf("      Work function        %g eV\n", m_work);
            printf("      Cluster size         %d\n", m_nsize);
        }

        // Initial situation: starting position
        double x = x0;
        double y = y0;
        double z = z0;

        // Store the energy [MeV].
        double e = ekin0;
        // Total distance covered
        double dsum = 0.0;
        // Pool of unused energy
        double epool = 0.0;

        auto cols = m_generator->Generate(GetKineticEnergy(),
                                          x0, y0, z0,
                                          xdir, ydir, zdir);

        bool bNClustersReachLimit = false;

        for (auto &col : cols)
        {

            // Cluster generated by recoil ion
            auto pos_col = col.GetPosition();
            const double x_col = pos_col.X();
            const double y_col = pos_col.Y();
            const double z_col = pos_col.Z();
            const double ene_incident = col.GetIncidentEnergy();
            const double ene_recoil = col.GetRecoilEnergy();
            if (ene_recoil > 0 && IsInside(x_col, y_col, z_col))
            {
                cluster cluster_recoil;
                cluster_recoil.x = x_col;
                cluster_recoil.y = y_col;
                cluster_recoil.z = z_col;
                cluster_recoil.t = t0;

                cluster_recoil.electrons = std::round(ene_recoil / m_work);
                cluster_recoil.ec = ene_recoil;
                cluster_recoil.kinetic = ene_incident;

                if (NewClusterPushable())
                {

                    if (cluster_recoil.electrons > 0)
                    {
                        if (m_debug)
                        {
                            std::cout << hdr << "Cluster " << m_clusters.size() << "\n    at ("
                                      << cluster_recoil.x << ", " << cluster_recoil.y << ", " << cluster_recoil.z
                                      << "),\n    e = " << cluster_recoil.ec << ",\n    n = "
                                      << cluster_recoil.electrons << ",\n    pool = "
                                      << cluster_recoil.kinetic << " eV.\n";
                        }

                        m_clusters.push_back(cluster_recoil);
                    }
                }
                else
                {
                    bNClustersReachLimit = true;
                    break;
                }
            }

            // Clusters generated during this step
            double ene_step = ene_incident - ene_recoil;

            const double dr = col.GetDistanceToNextCollision();
            auto dir_col = col.GetScatteringDirection();
            const double dir_x_col = dir_col.X();
            const double dir_y_col = dir_col.Y();
            const double dir_z_col = dir_col.Z();

            const int nElectrons = std::round(col.GetEnergyLoss() / m_work);
            const int nClusters = m_nsize < 0 ? nElectrons : std::ceil(nElectrons / m_nsize);
            if (nClusters == 0 || nElectrons == 0)
                continue;

            const int nElectronsInCluster = nElectrons / nClusters;
            const double eneCluster = std::round(col.GetEnergyLoss() / nClusters);

            const int nDiv = nClusters + 1;

            for (int iCluster = 0; iCluster < nClusters; ++iCluster)
            {

                cluster newcluster;
                double x_cls, y_cls, z_cls;
                x_cls = x_col + dr * dir_x_col * (iCluster + 1) / nDiv;
                y_cls = y_col + dr * dir_y_col * (iCluster + 1) / nDiv;
                z_cls = z_col + dr * dir_z_col * (iCluster + 1) / nDiv;

                if (!IsInside(x_cls, y_cls, z_cls))
                    continue;

                newcluster.x = x_cls;
                newcluster.y = y_cls;
                newcluster.z = z_cls;
                newcluster.t = t0;

                newcluster.electrons = nElectronsInCluster;
                newcluster.ec = eneCluster;
                newcluster.kinetic = ene_step;

                ene_step -= eneCluster;

                if (NewClusterPushable())
                {
                    if (m_debug)
                    {
                        std::cout << hdr << "Cluster " << m_clusters.size() << "\n    at ("
                                  << newcluster.x << ", " << newcluster.y << ", " << newcluster.z
                                  << "),\n    e = " << newcluster.ec << ",\n    n = "
                                  << newcluster.electrons << ",\n    pool = "
                                  << newcluster.kinetic << " eV.\n";
                    }

                    m_clusters.push_back(newcluster);
                }
                else
                {
                    bNClustersReachLimit = true;
                    break;
                }
            }
        }

        if (bNClustersReachLimit)
        {
            std::cerr << hdr << "Exceeded maximum number of clusters.\n";
        }

        return true;
        // finished generating
    }

    bool TrackTrimSQLite::GetCluster(double &xcls, double &ycls, double &zcls,
                                     double &tcls, int &n, double &e, double &extra)
    {

        if (m_debug)
        {
            printf("Current cluster: %d, array size: %ld",
                   m_currcluster, m_clusters.size());
        }
        // Stop if we have exhausted the list of clusters.
        if (m_currcluster >= m_clusters.size())
            return false;

        xcls = m_clusters[m_currcluster].x;
        ycls = m_clusters[m_currcluster].y;
        zcls = m_clusters[m_currcluster].z;
        tcls = m_clusters[m_currcluster].t;

        n = m_clusters[m_currcluster].electrons;
        e = m_clusters[m_currcluster].ec;
        extra = m_clusters[m_currcluster].kinetic;
        // Move to next cluster
        ++m_currcluster;
        return true;
    }

    bool TrackTrimSQLite::IsInside(const double x, const double y, const double z)
    {
        Medium *medium = NULL;

        // Check that the cluster is in an ionisable medium and within bounding box
        if (!m_sensor->GetMedium(x, y, z, medium))
        {
            if (m_debug)
            {
                std::cout << "No medium at position ("
                          << x << "," << y << "," << z << ").\n";
            }
            return false;
        }
        else if (!medium->IsIonisable())
        {
            if (m_debug)
            {
                std::cout << "Medium at ("
                          << x << "," << y << "," << z << ") is not ionisable.\n";
            }
            return false;
        }
        else if (!m_sensor->IsInArea(x, y, z))
        {
            if (m_debug)
            {
                std::cout << "Cluster at ("
                          << x << "," << y << "," << z << ") outside bounding box.\n";
            }
            return false;
        }
        return true;
    }

}
