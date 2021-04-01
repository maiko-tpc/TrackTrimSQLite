#include "TrackGenerator.hpp"

#include <algorithm>

void TrackGeneratorA::EraseHighEnergyCollisionsFromTrack(CollisionCollection &_track,
                                                         double _ekin)
{
    auto it_ekin = std::find_if(_track.begin(), _track.end(),
                                [&_ekin](Collision &_col) -> bool {
                                    double e = _col.GetIncidentEnergy();
                                    double e_next = e - _col.GetEnergyLoss() - _col.GetRecoilEnergy();
                                    return (e_next <= _ekin && _ekin <= e);
                                });

    // Incident energy at 0 th collision is smaller than input _ekin.
    if (it_ekin == _track.end())
        throw std::runtime_error("Kinetic energy out of range.");

    // Remove collisions before ekin
    if (it_ekin != _track.begin())
        _track.erase(_track.begin(), it_ekin);
}

TrackGenerator::Transform TrackGeneratorA::SetFirstCollision(CollisionIterator _it,
                                                             double _ekin, double _x, double _y, double _z,
                                                             double _dx, double _dy, double _dz) const
{

    Transform mat;
    //kinetic energy after recoil
    double e = _it->GetIncidentEnergy() - _it->GetRecoilEnergy();
    double de = _it->GetEnergyLoss();
    double dr = _it->GetDistanceToNextCollision();

    //update distance to next collision
    double de_ratio = de != 0 ? (de - (e - _ekin)) / de : 1;
    double dr_update = dr * de_ratio;
    double de_update = de * de_ratio;
    _it->SetIncidentEnergy(_ekin);
    _it->SetRecoilEnergy(0); //<-- initial vector
    _it->SetEnergyLoss(de_update);
    _it->SetDistanceToNextCollision(dr_update);

    //rotate direction
    Vector dx0 = _it->GetIncidentDirection();
    Mm::Normalize(dx0);
    Vector dx1 = _it->GetScatteringDirection();
    Mm::Normalize(dx1);
    Vector dx_inc(_dx, _dy, _dz);
    Mm::Normalize(dx_inc);
    //mat = MakeRotaionMatrix(axis, Angle(dx0, dx_inc));
    if (Mm::Norm(dx1) != 0)
    {
        mat = Mm::MakeRotaionMatrix(dx1, dx_inc);
    }
    else if (Mm::Norm(dx0) != 0)
    {
        mat = Mm::MakeRotaionMatrix(dx0, dx_inc);
    }
    else
    {
        // Assume track sample at 0th collision is parallel to x-axis
        Vector dx_initial(1, 0, 0);
        mat = Mm::MakeRotaionMatrix(dx_initial, dx_inc);
    }
    //update
    _it->SetPosition(_x, _y, _z);
    dx0 = dx_inc;
    _it->SetIncidentDirection(dx0);

    //dx1 = mat.Apply(dx1);
    //Mm::Normalize(dx1);
    if (Mm::Norm(dx1) != 0)
        dx1 = dx_inc;
    _it->SetScatteringDirection(dx1);
    return mat;
}

// Adapt collision at _it to the collision just before
// _it must not be the first iterator
void TrackGeneratorA::AdaptCollision(CollisionIterator _it, Transform _mat)
{
    auto it_before = _it - 1;
    auto pos = Mm::Add(it_before->GetPosition(),
                       Mm::Scale(it_before->GetDistanceToNextCollision(),
                                 it_before->GetScatteringDirection()));
    auto dx1_before = it_before->GetScatteringDirection();

    auto dx0 = _it->GetIncidentDirection();

    auto dx1 = _it->GetScatteringDirection();

    dx0 = _mat.Apply(dx0);
    Mm::Normalize(dx0);
    //std::cout << "conpare -> " << dx1_before.Y() << " vs " << dx0.Y() << " " << std::endl;
    dx1 = _mat.Apply(dx1);
    Mm::Normalize(dx1);
    _it->SetPosition(pos);
    _it->SetIncidentDirection(dx0);
    _it->SetScatteringDirection(dx1);
};

void TrackGeneratorB::PhiRotationRandom(CollisionIterator _it, Transform &_mat)
{
    auto dx0 = _it->GetIncidentDirection();
    Mm::Normalize(dx0);
    auto dx1 = _it->GetScatteringDirection();
    Mm::Normalize(dx1);
    double dp_before = Mm::DotProduct(dx0, dx1);
    if (Mm::Norm(dx0) != 0)
    {
        auto angle_rand = 2 * M_PI * Random();
        //auto angle_rand = 1;
        auto rotation_rand = Mm::MakeRotaionMatrix(dx0, angle_rand);
        //Vector x(1, 0, 0);
        //auto rotation_rand = MakeRotaionMatrix(x, angle_rand);
        dx1 = rotation_rand.Apply(dx1);
        Mm::Normalize(dx1);
        _it->SetScatteringDirection(dx1);
        _mat = rotation_rand * _mat;
        _mat.Normalize();
        double dp_after = Mm::DotProduct(dx0, dx1);
        //std::cout << dp_before << " -> " << dp_after << std::endl;
    }
};

TrackGenerator::CollisionCollection TrackGeneratorC::GetTransferDestinationCandidates(double _ene,
                                                                                      double _ene_min, double _ene_max)
{
    CollisionDBHandler db(GetFileName());
    int nTracks = db.GetNumberOfTracks();
    std::string sConstraint;
    // 1 : Energy after collision is nearly-equal to this collision
    sConstraint += std::to_string(_ene_min) + " <= e_inc - e_rec AND ";
    sConstraint += "e_inc - e_rec  <= " + std::to_string(_ene_max) + " AND ";
    // 2 : Energy at next Collision is smaller than this collision
    sConstraint += "e_inc - e_rec - de < " + std::to_string(_ene) + " AND ";
    // 3 : Not the last collision (== energy after collision is NOT 0)
    sConstraint += "0 < e_inc - e_rec - de AND ";
    // 4 : Scattering direction is defined
    sConstraint += "0 < dr ";
    if (GetTrackIDMin() >= 0)
        sConstraint += " AND " + std::to_string(GetTrackIDMin()) + " <= track_id ";
    if (GetTrackIDMax() >= 0)
        sConstraint += " AND track_id <= " + std::to_string(GetTrackIDMax());

    return db.GetCollisions(sConstraint);
}

bool TrackGeneratorC::Transfer(CollisionCollection &_track, CollisionIterator &_it, Transform &_mat)
{
    // Kinetic energy after current collision
    double ene = _it->GetIncidentEnergy() - _it->GetRecoilEnergy();
    // Energy loss before next collision in current collision sequence
    double de = _it->GetEnergyLoss();
    // Distance to next collision in current collision sequence
    double dr = _it->GetDistanceToNextCollision();

    // Last collision
    if (ene - de == 0)
    {
        std::cout << " Attempt to transfer at end of the track -> ignore" << std::endl;
        return false;
    }

    // No energy loss
    if (de == 0)
    {
        std::cout << " No energy loss before next collision -> ignore" << std::endl;
        return false;
    }

    double ene_min = ene - de * GetEnergyMarginRatio();
    double ene_max = ene + de * GetEnergyMarginRatio();

    // Candidate : Collisions whose kinetic energy are nearly equal to that of current one.
    //   -> Collision next to the candidate is connected to the current collision
    CollisionCollection cols = GetTransferDestinationCandidates(ene, ene_min, ene_max);

    // Never be called ?
    if (cols.size() == 0)
    {
        std::cerr << "TrackGenerator::No collision for transfer found" << std::endl;
        //std::cerr << "energy left " << ene - de << std::endl;
        //throw std::runtime_error("No collision for transfer found");
        return false;
    }

    // Randomly select one collision to transfer
    auto col_selected = cols.begin() + GetRandomInteger(0, cols.size() - 1);
    int trackID_selected = col_selected->GetTrackID();
    int collisionID_selected = col_selected->GetCollisionID();

    // Incident energy of destination collision
    double e_inc_transfer = col_selected->GetIncidentEnergy();
    e_inc_transfer -= col_selected->GetRecoilEnergy();
    e_inc_transfer -= col_selected->GetEnergyLoss();

    // Retrieve collision sequence to be connected
    CollisionCollection cols_transfer = GetTrack(trackID_selected);
    auto it_begin_transfer = std::find_if(cols_transfer.begin(),
                                          cols_transfer.end(),
                                          [&collisionID_selected](Collision &_col) -> bool {
                                              return _col.GetCollisionID() == collisionID_selected;
                                          });
    //Never be called
    if (it_begin_transfer == cols_transfer.end())
    {
        throw std::runtime_error("Collision to transfer is not recorded ???");
    }
    // Collision similar to current one
    it_begin_transfer = cols_transfer.erase(cols_transfer.begin(), it_begin_transfer);
    // Collision to be connected
    ++it_begin_transfer;

    // DO UPDATE
    // Update current collision
    double de_ratio = de != 0 ? (ene - e_inc_transfer) / de : 1;
    double dr_update = dr * de_ratio;
    double de_update = de * de_ratio;

    std::cout << "    Transfer " << de_ratio << std::endl;
    std::cout << "             " << ene << " " << de << " " << dr << std::endl;
    std::cout << "             " << e_inc_transfer << " " << de_update << " " << dr_update << std::endl;

    double de_update_new = ene - e_inc_transfer;

    double dr_update_new;
    double dr_ratio_new;
    {
        double de_src = de; // non-zero
        double de_dst = cols_transfer.begin()->GetEnergyLoss();

        double det = std::abs(de_src - de_update_new) - std::abs(de_dst - de_update_new);

        auto distance = [](double a, double b) -> double { return std::abs(std::log(a) - std::log(b)); };
        // de_dst_is closer to de_update_new
        //       -> dr is calculated based on dr/dE of transfer destination.
        std::cout << de_update_new << " <- which is closer ? " << de_src << " " << de_dst << std::endl;
        if (de_dst > 0 && distance(de_dst, de_update_new) < distance(de_src, de_update_new))
        {
            double dr_dst = cols_transfer.begin()->GetDistanceToNextCollision();
            dr_ratio_new = dr_dst / de_dst;
        }
        // de_src_is closer to de_update_new or de_dst == 0
        //       -> dr is calculated based on dr / dE of transfer source.
        else
        {
            double dr_src = dr;
            dr_ratio_new = dr_src / de_src;
        }
    }
    dr_update_new = de_update_new * dr_ratio_new;

    std::cout << "             " << dr_update << " OR " << dr_update_new << " ratio " << dr_ratio_new << std::endl;
    //std::cout << "             " << de_update << " OR " << de_update_new << std::endl;

    //_it->SetEnergyLoss(de_update);
    //_it->SetDistanceToNextCollision(dr_update);

    _it->SetEnergyLoss(de_update_new);
    _it->SetDistanceToNextCollision(dr_update_new);

    // Save position of current collision to restore _it after data migration
    auto current_position = _it - _track.begin();

    // Erase collisions after current collision (_it is destroyed.)
    _track.erase(_it + 1, _track.end());
    // Number of collisions to be add
    auto n_collisions_after_transfer = cols_transfer.end() - it_begin_transfer;
    // Reserve data buffer to store collisions transferred to
    _track.reserve(_track.size() + n_collisions_after_transfer);
    // Copy collisions
    std::copy(it_begin_transfer, cols_transfer.end(), std::back_inserter(_track));
    // Reallocate current collision to _it
    _it = _track.begin() + current_position;

    // Update rotation matrix
    auto dir_incident_transfer = (_it + 1)->GetIncidentDirection();
    auto dir_scattering = _it->GetScatteringDirection();
    // std::cout << "    norm = " << Norm(dir_incident_transfer) << " " << Norm(dir_scattering) << std::endl;
    //std::cout << "A " << Mm::Norm(dir_incident_transfer) << " " << Mm::Norm(dir_scattering) << std::endl;
    auto mat_transfer = Mm::MakeRotaionMatrix(dir_incident_transfer, dir_scattering);
    //std::cout << "B" << std::endl;
    _mat = mat_transfer;

    return true;
}