#include "SrimDataFile.hpp"
#include "FundamentalConstants.hpp"


namespace Mm{
SrimDataFile::SrimDataFile(){};
SrimDataFile::~SrimDataFile(){};

bool SrimDataFile::ReadFile(const std::string &_filename)
{
    // SRMREA
    std::string fClassName("SrimDataFile");
  const std::string hdr = fClassName + "::ReadFile:\n    ";
  // Open the material list.
  std::ifstream fsrim;
  fsrim.open(_filename.c_str(), std::ios::in);
  if (fsrim.fail()) {
    std::cerr << hdr << "Could not open SRIM  file " << _filename
              << " for reading.\n    The file perhaps does not exist.\n";
    return false;
  }
  unsigned int nread = 0;

  // Read the header
  if (fDebug) {
    std::cout << hdr << "SRIM header records from file " << _filename << "\n";
  }
  const int size = 100;
  char line[size];
  while (fsrim.getline(line, 100, '\n')) {
    nread++;
    if (strstr(line, "SRIM version") != NULL) {
      if (fDebug) std::cout << "\t" << line << "\n";
    } else if (strstr(line, "Calc. date") != NULL) {
      if (fDebug) std::cout << "\t" << line << "\n";
    } else if (strstr(line, "Ion =") != NULL) {
      break;
    }
  }

  // Identify the ion
  char* token = NULL;
  token = strtok(line, " []=");
  token = strtok(NULL, " []=");
  token = strtok(NULL, " []=");
  // Set the ion charge.
  fZ = std::atof(token);
  //m_chargeset = true;
  token = strtok(NULL, " []=");
  token = strtok(NULL, " []=");
  token = strtok(NULL, " []=");
  // Set the ion mass (convert amu to eV).
  fMass = std::atof(token) * AtomicMassUnitElectronVolt;

  // Find the target density
  if (!fsrim.getline(line, 100, '\n')) {
    std::cerr << hdr << "Premature EOF looking for target density (line "
              << nread << ").\n";
    return false;
  }
  nread++;
  if (!fsrim.getline(line, 100, '\n')) {
    std::cerr << hdr << "Premature EOF looking for target density (line "
              << nread << ").\n";
    return false;
  }
  nread++;
  token = strtok(line, " ");
  token = strtok(NULL, " ");
  token = strtok(NULL, " ");
  token = strtok(NULL, " ");
  //SetDensity(std::atof(token));
  fDensity = std::atof(token);

  // Check the stopping units
  while (fsrim.getline(line, 100, '\n')) {
    nread++;
    if (strstr(line, "Stopping Units") == NULL) continue;
    if (strstr(line, "Stopping Units =  MeV / (mg/cm2)") != NULL) {
      if (fDebug) {
        std::cout << hdr << "Stopping units: MeV / (mg/cm2) as expected.\n";
      }
      break;
    }
    std::cerr << hdr << "Unknown stopping units. Aborting (line " << nread
              << ").\n";
    return false;
  }

  // Skip to the table
  while (fsrim.getline(line, 100, '\n')) {
    nread++;
    if (strstr(line, "-----------") != NULL) break;
  }

  // Read the table line by line
  fvEkin.clear();
  fvEmLoss.clear();
  fvHdLoss.clear();
  fvRange.clear();
  fvTransverseStraggling.clear();
  fvLongitudinalStraggling.clear();
  unsigned int ntable = 0;
  while (fsrim.getline(line, 100, '\n')) {
    nread++;
    if (strstr(line, "-----------") != NULL) break;
    // Energy
    token = strtok(line, " ");
    fvEkin.push_back(atof(token));
    token = strtok(NULL, " ");
    if (strcmp(token, "eV") == 0) {
      fvEkin[ntable] *= 1.0e-6;
    } else if (strcmp(token, "keV") == 0) {
      fvEkin[ntable] *= 1.0e-3;
    } else if (strcmp(token, "GeV") == 0) {
      fvEkin[ntable] *= 1.0e3;
    } else if (strcmp(token, "MeV") != 0) {
      std::cerr << hdr << "Unknown energy unit " << token << "; aborting\n";
      return false;
    }
    // EM loss
    token = strtok(NULL, " ");
    fvEmLoss.push_back(atof(token));
    // HD loss
    token = strtok(NULL, " ");
    fvHdLoss.push_back(atof(token));
    // Projected range
    token = strtok(NULL, " ");
    fvRange.push_back(atof(token));
    token = strtok(NULL, " ");
    if (strcmp(token, "A") == 0) {
      fvRange[ntable] *= 1.0e-8;
    } else if (strcmp(token, "um") == 0) {
      fvRange[ntable] *= 1.0e-4;
    } else if (strcmp(token, "mm") == 0) {
      fvRange[ntable] *= 1.0e-1;
    } else if (strcmp(token, "m") == 0) {
      fvRange[ntable] *= 1.0e2;
    } else if (strcmp(token, "km") == 0) {
      fvRange[ntable] *= 1.0e5;
    } else if (strcmp(token, "cm") != 0) {
      std::cerr << hdr << "Unknown distance unit " << token << "; aborting\n";
      return false;
    }
    // Longitudinal straggling
    token = strtok(NULL, " ");
    fvLongitudinalStraggling.push_back(atof(token));
    token = strtok(NULL, " ");
    if (strcmp(token, "A") == 0) {
      fvLongitudinalStraggling[ntable] *= 1.0e-8;
    } else if (strcmp(token, "um") == 0) {
      fvLongitudinalStraggling[ntable] *= 1.0e-4;
    } else if (strcmp(token, "mm") == 0) {
      fvLongitudinalStraggling[ntable] *= 1.0e-1;
    } else if (strcmp(token, "m") == 0) {
      fvLongitudinalStraggling[ntable] *= 1.0e2;
    } else if (strcmp(token, "km") == 0) {
      fvLongitudinalStraggling[ntable] *= 1.0e5;
    } else if (strcmp(token, "cm") != 0) {
      std::cerr << hdr << "Unknown distance unit " << token << "; aborting\n";
      return false;
    }
    // Transverse straggling
    token = strtok(NULL, " ");
    fvTransverseStraggling.push_back(atof(token));
    token = strtok(NULL, " ");
    if (strcmp(token, "A") == 0) {
      fvTransverseStraggling[ntable] *= 1.0e-8;
    } else if (strcmp(token, "um") == 0) {
      fvTransverseStraggling[ntable] *= 1.0e-4;
    } else if (strcmp(token, "mm") == 0) {
      fvTransverseStraggling[ntable] *= 1.0e-1;
    } else if (strcmp(token, "m") == 0) {
      fvTransverseStraggling[ntable] *= 1.0e2;
    } else if (strcmp(token, "km") == 0) {
      fvTransverseStraggling[ntable] *= 1.0e5;
    } else if (strcmp(token, "cm") != 0) {
      std::cerr << hdr << "Unknown distance unit " << token << "; aborting\n";
      return false;
    }

    // Increment table line counter
    ++ntable;
  }

  // Find the scaling factor and convert to MeV/cm
  double scale = -1.;
  while (fsrim.getline(line, 100, '\n')) {
    nread++;
    if (strstr(line, "=============") != NULL) {
      break;
    } else if (strstr(line, "MeV / (mg/cm2)") != NULL ||
               strstr(line, "MeV/(mg/cm2)") != NULL) {
      token = strtok(line, " ");
      scale = std::atof(token);
    }
  }
  if (scale < 0) {
    std::cerr << hdr << "Did not find stopping unit scaling; aborting.\n";
    return false;
  }
  scale *= 1.e3;
  for (unsigned int i = 0; i < ntable; ++i) {
    fvEmLoss[i] *= scale;
    fvHdLoss[i] *= scale;
  }

  // Seems to have worked
  if (fDebug) {
    std::cout << hdr << "Successfully read " << _filename << "(" << nread
              << " lines).\n";
  }
  return true;
}

void SrimDataFile::Show()const
{
  for(int i = 0; i< fvEkin.size(); ++ i)
  {
    std::cout << fvEkin.at(i) << " "
              << fvEmLoss.at(i) << " "
              << fvHdLoss.at(i) << " "
              << fvRange.at(i) << " "
              << fvLongitudinalStraggling.at(i) <<  " "
              << fvTransverseStraggling.at(i) << std::endl;
  }
}

}