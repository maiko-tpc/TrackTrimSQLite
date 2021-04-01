#pragma once

namespace Mm
{
    template <typename T>
    struct Matrix
    {
        Matrix()
            : fXX(1), fXY(0), fXZ(0),
              fYX(0), fYY(1), fYZ(0),
              fZX(0), fZY(0), fZZ(1){};

        double fXX, fXY, fXZ;
        double fYX, fYY, fYZ;
        double fZX, fZY, fZZ;

        T Apply(const T &_in)
        {
            double x = fXX * _in.X() + fXY * _in.Y() + fXZ * _in.Z();
            double y = fYX * _in.X() + fYY * _in.Y() + fYZ * _in.Z();
            double z = fZX * _in.X() + fZY * _in.Y() + fZZ * _in.Z();
            return T(x, y, z);
        }

        T Row(int _i) const
        {
            if (_i == 0)
                return T(fXX, fXY, fXZ);
            else if (_i == 1)
                return T(fYX, fYY, fYZ);
            else if (_i == 2)
                return T(fZX, fZY, fZZ);
            else
                throw std::out_of_range("Matrix::Row");
        }

        T Column(int _i) const
        {
            if (_i == 0)
                return T(fXX, fYX, fZX);
            else if (_i == 1)
                return T(fXY, fYY, fZY);
            else if (_i == 2)
                return T(fXZ, fYZ, fZZ);
            else
                throw std::out_of_range("Matrix::Column");
        }

        void Initialize()
        {
            fXX = fYY = fZZ = 1;
            fXY = fXZ = 0;
            fYX = fYZ = 0;
            fZX = fZY = 0;
        }

        double Determinant() const
        {
            double pos = fXX * fYY * fZZ + fYX * fZY * fXZ + fZX * fXY * fYZ;
            double neg = fXZ * fYY * fZX + fYZ * fZY * fXX + fZZ * fXY * fYX;
            return pos - neg;
        }

        void Scale(double _sca)
        {
            fXX = _sca * fXX;
            fXY = _sca * fXY;
            fXZ = _sca * fXZ;

            fYX = _sca * fYX;
            fYY = _sca * fYY;
            fYZ = _sca * fYZ;

            fZX = _sca * fZX;
            fZY = _sca * fZY;
            fZZ = _sca * fZZ;
        }

        void Normalize()
        {
            auto det = Determinant();
            if (det != 0)
            {
                Scale(1 / std::cbrt(det));
            }
        }
    };

    template <typename T>
    T Scale(double _sc, const T &_vec)
    {
        return T(_sc * _vec.X(), _sc * _vec.Y(), _sc * _vec.Z());
    }

    template <typename T>
    T Add(const T &_l, const T &_r)
    {
        return T(_l.X() + _r.X(), _l.Y() + _r.Y(), _l.Z() + _r.Z());
    }

    template <typename T>
    double DotProduct(const T &_l, const T &_r)
    {
        return _l.X() * _r.X() + _l.Y() * _r.Y() + _l.Z() * _r.Z();
    }

    template <typename T>
    double Norm(const T &_v)
    {
        return std::sqrt(DotProduct(_v, _v));
    }

    template <typename T>
    double Angle(const T &_l, const T &_r)
    {
        auto cross = DotProduct(_l, _r);
        double normalization = Norm(_l) * Norm(_r);
        if (normalization == 0)
            throw std::runtime_error("Angle(): zero division.");
        //cos (theta)
        double cs = cross / normalization;
        return std::acos(cs);
    }

    template <typename T>
    bool Normalize(T &_v)
    {
        double norm = Norm(_v);
        if (norm == 0)
            return false;
        _v = T(_v.X() / norm, _v.Y() / norm, _v.Z() / norm);
        return true;
    }

    template <typename T>
    T CrossProduct(const T &_l, const T &_r)
    {
        double x = _l.Y() * _r.Z() - _l.Z() * _r.Y();
        double y = _l.Z() * _r.X() - _l.X() * _r.Z();
        double z = _l.X() * _r.Y() - _l.Y() * _r.X();

        return T(x, y, z);
    }

    template <typename T>
    Matrix<T> operator*(const Matrix<T> &_l, const Matrix<T> &_r)
    {
        Matrix<T> ret;

        ret.fXX = DotProduct(_l.Row(0), _r.Column(0));
        ret.fXY = DotProduct(_l.Row(0), _r.Column(1));
        ret.fXZ = DotProduct(_l.Row(0), _r.Column(2));

        ret.fYX = DotProduct(_l.Row(1), _r.Column(0));
        ret.fYY = DotProduct(_l.Row(1), _r.Column(1));
        ret.fYZ = DotProduct(_l.Row(1), _r.Column(2));

        ret.fZX = DotProduct(_l.Row(2), _r.Column(0));
        ret.fZY = DotProduct(_l.Row(2), _r.Column(1));
        ret.fZZ = DotProduct(_l.Row(2), _r.Column(2));

        return ret;
    }

    template <typename T>
    Matrix<T> MakeRotaionMatrix(const T &_axis, double _degree)
    {

        double norm = Norm(_axis);
        Matrix<T> mat;
        if (norm == 0)
        {
            throw std::runtime_error("MakeRotationMatrix() : Zero vector input !");
        }

        double nx = _axis.X() / norm;
        double ny = _axis.Y() / norm;
        double nz = _axis.Z() / norm;
        double sn = sin(_degree);
        double cs = cos(_degree);

        mat.fXX = nx * nx * (1 - cs) + cs;
        mat.fXY = nx * ny * (1 - cs) - nz * sn;
        mat.fXZ = nz * nx * (1 - cs) + ny * sn;

        mat.fYX = nx * ny * (1 - cs) + nz * sn;
        mat.fYY = ny * ny * (1 - cs) + cs;
        mat.fYZ = ny * nz * (1 - cs) - nx * sn;

        mat.fZX = nz * nx * (1 - cs) - ny * sn;
        mat.fZY = ny * nz * (1 - cs) + nx * sn;
        mat.fZZ = nz * nz * (1 - cs) + cs;

        return mat;
    };

    //Make rotation matrix which transform _dir1 to _dir2
    template <typename T>
    Matrix<T> MakeRotaionMatrix(const T &_dir1, const T &_dir2)
    {
        T dir1 = _dir1;
        Normalize(dir1);
        T dir2 = _dir2;
        Normalize(dir2);
        auto axis = CrossProduct(dir1, dir2);
        Normalize(axis);
        return MakeRotaionMatrix(axis, Angle(dir1, dir2));
    };
}
