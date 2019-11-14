/*==============================================================================

                                    O  F  E  L  I

                            Object  Finite  Element  Library

  ==============================================================================

   Copyright (C) 1998 - 2020 Rachid Touzani

   This file is part of OFELI.

   OFELI is free software: you can redistribute it and/or modify
   it under the terms of the GNU Lesser General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   OFELI is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU Lesser General Public License for more details.

   You should have received a copy of the GNU Lesser General Public License
   along with OFELI. If not, see <http://www.gnu.org/licenses/>.

  ==============================================================================

             Definition of class 'IPF' to manage Input Project files

  ==============================================================================*/

#ifndef __IPF_H
#define __IPF_H

#include "OFELI_Config.h"
#include <complex>
using std::complex;

namespace OFELI {
/*!
 *  \addtogroup OFELI
 *  @{
 */

/*! \file IPF.h
 *  \brief Definition file for class IPF.
 */

/** \defgroup Util Utilities
 *  \brief Utility functions and classes
 */


/** \class IPF
 * \ingroup IO
 *  \brief To read project parameters from a file in IPF format.
 *
 * \details This class can be used to acquire various parameters from a
 * parameter file of IPF (Input Project File).
 * The declaration of an instance of this class avoids reading data
 * in your main program. The acquired parameters are retrieved
 * through information members of the class. Note that all the
 * parameters have default values
 *
 * \author Rachid Touzani
 * \copyright GNU Lesser Public License
 */

class IPF {

 public:

#ifndef DOXYGEN_SHOULD_SKIP_THIS
   struct mat_prop {
       string _density, _electric_cond, _electric_perm, _magnetic_perm;
       string _poisson, _thermal_cond, _rho_cp, _visc, _young;
   };
#endif /* DOXYGEN_SHOULD_SKIP_THIS */

/// \brief Default constructor
    IPF();

/// \brief Constructor that gives the data file name.
/// \details It reads parameters in IPF Format from this file.
    IPF(const string& file);

/** \brief Constructor that reads parameters in file <tt>file</tt> and
 *  prints header information for the calling program <tt>prog</tt>.
 *  It reads parameters in <tt>IPF</tt> Format from this file.
 */
    IPF(const string& prog,
        const string& file);

/// \brief Destructor
    ~IPF();

/// \brief Display acquired parameters
    real_t getDisplay();

/// \brief Return parameter read using keyword \b Verbose
    int getVerbose() { return _verbose; }

/// \brief Return parameter read using keyword \b Output.
/// \details This parameter can be used to control output behavior in a program.
    int getOutput() const { return _output; }

/// \brief Return parameter read using keyword \b Save.
/// \details This parameter can be used to control result saving in a program
/// (\e e.g. for a restarting purpose).
    int getSave() const { return _save; }

/// \brief Return parameter read using keyword \b Plot.
/// \details This parameter can be used to control result saving for plotting in a program.
    int getPlot() const { return _plot; }

/// \brief Return parameter read using keyword \b BC.
/// \details This parameter can be used to set a boundary condition flag.
    int getBC() const { return _bc; }

/// \brief Return parameter read using keyword \b BF.
/// \details This parameter can be used to set a body force flag.
    int getBF() const { return _bf; }

/// \brief Return parameter read using keyword \b SF.
/// \details This parameter can be used to set a surface force flag.
    int getSF() const { return _sf; }

/// \brief Return parameter read using keyword \b Init.
/// \details This parameter can be used to set an initial data flag.
    int getInit() const { return _ini; }

/// \brief Return parameter read using keyword \b Data.
/// \details This parameter can be used to set a various data flag.
    int getData() const { return _data; }

/// \brief Return parameter read using keyword \b NbSteps.
/// \details This parameter can be used to read a number of time steps.
    size_t getNbSteps() const { return _nb_steps; }

/// \brief Return parameter read using keyword \b NbIter.
/// \details This parameter can be used to read a number of iterations.
    size_t getNbIter() const { return _nb_iter; }

/// \brief Return parameter read using keyword \b TimeStep.
/// \details This parameter can be used to read a time step value.
    real_t getTimeStep() const { return _time_step; }

/// \brief Return parameter read using keyword \b MaxTime.
/// \details This parameter can be used to read a maximum time value.
    real_t getMaxTime() const { return _max_time; }

/// \brief Return parameter read using keyword \b Tolerance.
/// \details This parameter can be used to read a tolerance value to control convergence.
    real_t getTolerance() const { return _tolerance; }

/// \brief Return <tt>n</tt>-th parameter read using keyword <tt>IntPar</tt>
/// \details Here we have at most 20 integer extra parameters that can be used for any purpose.
/// Default value for <tt>n</tt> is <tt>1</tt>
    int getIntPar(size_t n=1) const { return _int_par[n-1]; }

/// \brief Return \a n-th parameter read using keyword \b StringPar.
/// \details Here we have at most 20 integer extra parameters that can be used for any purpose.
/// Default value for <tt>n</tt> is <tt>1</tt>
    string getStringPar(size_t n=1) const { return _string_par[n-1]; }

/// \brief Return <tt>n</tt>-th parameter read using keyword \b DoublePar.
/// \details Here we have at most 20 integer extra parameters that can be used for any purpose.
/// Default value for <tt>n</tt> is <tt>1</tt>
    real_t getDoublePar(size_t n=1) const { return _real_par[n-1]; }

/// \brief Return <tt>n</tt>-th parameter read using keyword \b PointDoublePar.
/// \details Here we have at most 20 integer extra parameters that can be used for any purpose.
/// Default value for <tt>n</tt> is <tt>1</tt>
    Point<real_t> getPointDoublePar(size_t n=1) const { return _point_double_par[n-1]; }

/// \brief Return <tt>n</tt>-th parameter read using keyword \b StringPar.
/// \details Here we have at most 20 integer extra parameters that can be used for any purpose.
/// Default value for <tt>n</tt> is <tt>1</tt>
    complex_t getComplexPar(size_t n=1) const { return _complex_par[n-1]; }

/// \brief Return parameter corresponding to a given label, when its value is a string
/// @param [in] label Label that identifies the string (read from input file)
/// If this label is not found an error message is displayed and program stops
    string getString(const string& label) const;

/** \brief Return parameter corresponding to a given label, when its value is a string
 *  \details Case where a default value is provided
 *  @param [in] label Label that identifies the string (read from input file)
 *  @param [in] def Default value: Value to assign if the sought parameter is not found
 */
    string getString(const string& label,
                           string  def) const;

/// \brief Return parameter corresponding to a given label, when its value is an integer
/// @param [in] label Label that identifies the integer number (read from input file)
/// If this label is not found an error message is displayed and program stops
    int getInteger(const string& label) const;

/** \brief Return parameter corresponding to a given label, when its value is an integer.
 *  \details Case where a default value is provided
 *  @param [in] label Label that identifies the integer number (read from input file).
 *  @param [in] def Default value: Value to assign if the sought parameter is not found
 */
    int getInteger(const string& label,
                         int     def) const;

/// \brief Return parameter corresponding to a given label, when its value is a real_t
/// @param [in] label Label that identifies the real number (read from input file).
/// If this label is not found an error message is displayed and program stops.
    real_t getDouble(const string& label) const;
   
/** \brief Return parameter corresponding to a given label, when its value is a real_t
 *  \details Case where a default value is provided
 *  @param [in] label Label that identifies the real number (read from input file)
 *  @param [in] def Default value: Value to assign if the sought parameter is not found
 */
    real_t getDouble(const string& label,
                     real_t        def) const;

/// \brief Return parameter corresponding to a given label, when its value is a complex number
/// @param [in] label Label that identifies the complex number (read from input file)
/// If this label is not found an error message is displayed and program stops
    complex_t getComplex(const string& label) const;
   
/** \brief Return parameter corresponding to a given label, when its value is a complex number
 *  \details Case where a default value is provided
 *  @param [in] label Label that identifies the complex number (read from input file)
 *  @param [in] def Default value: Value to assign if the sought parameter is not found
 */
    complex_t getComplex(const string&   label,
                         complex_t       def) const;

/** \brief check if the project file contains a given parameter
 *  @param [in] label Label that identifies the label to seek in file
 *  @return <tt>0</tt> if the parameter is not found,
 *          <tt>n</tt> if the parameter is found, where <tt>n</tt> is the parameter index
 *          in the parameter list
 */
    int contains(const string& label) const;

/** \brief Read an array of real values, corresponding to a given label
 *  @param [in] label Label that identifies the array (read from input file).
 *  @param [in] a Vector that contain the array. The vector is properly resized before filling.
 *  @remark If this label is not found an error message is displayed.
 */
    void get(const string& label,
             Vect<real_t>& a) const;

/** \brief Return an array entry for a given label
 *  @param [in] label Label that identifies the array (read from input file).
 *  @param [in] j Index of entry in the array (Starting from 1)
 *  @remark If this label is not found an error message is displayed and program stops.
 */
    real_t getArraySize(const string& label,
                        size_t        j) const;

/** \brief Return integer parameter corresponding to a given label
 *  @param [in] label Label that identifies the integer number (read from input file).
 *  @param [out] a Returned value.
 *  If this label is not found an error message is displayed and program stops.
 *  Note: This member function can be used instead of getInteger
 */
    void get(const string& label,
             int&          a)    const;

/** \brief Return real parameter corresponding to a given label
 *  @param [in] label Label that identifies the real (real_t) number (read from input file).
 *  @param [out] a Returned value.
 *  If this label is not found an error message is displayed and program stops.
 *  Note: This member function can be used instead of getReal_T
 */
    void get(const string& label,
             real_t&       a) const;

/** \brief Return complex parameter corresponding to a given label
 *  @param [in] label Label that identifies the complex number (read from input file).
 *  @param [out] a Returned value.
 *  If this label is not found an error message is displayed and program stops.
 */
    void get(const string& label,
             complex_t&    a) const;

/** \brief Return string parameter corresponding to a given label
 *  @param [in] label Label that identifies the atring (read from input file).
 *  @param [out] a Returned value.
 *  Note: This member function can be used instead of getString
 *  If this label is not found an error message is displayed and program stops.
 *  Note: This member function can be used instead of getString
 */
    void get(const string& label,
             string&       a) const;

/// \brief Return parameter read using keyword \b Project.
/// \details This parameter can be used to read a project's name.
    string getProject() const { return _project; }
   
/// \brief Return pameter using keyword \b Mesh
    string getDomainFile() const { return _domain_file; }
   
/// \brief Return <tt>i</tt>-th parameter read using keyword \b mesh_file
/// \details Here we have at most 10 integer extra parameters that can be used for any purpose.
/// Default value for <tt>i</tt> is <tt>1</tt>
    string getMeshFile(size_t i=1) const { return _mesh_file[i-1]; }

/// \brief Return parameter read using keyword \b InitFile.
/// \details This parameter can be used to read an initial data file name.
    string getInitFile() const { return _init_file; }

/// \brief Return parameter read using keyword \b RestartFile.
/// \details This parameter can be used to read a restart file name.
    string getRestartFile() const { return _restart_file; }

/// \brief Return parameter read using keyword \b BCFile.
/// \details This parameter can be used to read a boundary condition file name.
    string getBCFile() const { return _bc_file; }

/// \brief Return parameter read using keyword \b BFFile.
/// \details This parameter can be used to read a body force file name.
    string getBFFile() const { return _bf_file; }

/// \brief Return parameter read using keyword \b SFFile.
/// \details This parameter can be used to read a source force file name.
    string getSFFile() const { return _sf_file; }

/// \brief Return parameter read using keyword \b SaveFile.
/// \details This parameter can be used to read a save file name.
    string getSaveFile() const { return _save_file; }

/// \brief Return <tt>i</tt>-th parameter read using keyword \b PlotFile.
/// \details Here we have at most 10 integer extra parameters that can be used for plot file names.
/// Default value for <tt>i</tt> is <tt>1</tt>
    string getPlotFile(int i=1) const { return _plot_file[i-1]; }

/// \brief Return parameter read using keyword \b DataFile.
/// \details This parameter can be used to read a Prescription file.
    string getPrescriptionFile(int i=1) const { return _data_file[i-1]; }
#ifndef DOXYGEN_SHOULD_SKIP_THIS
    string getDataFile(int i=1) const { return _data_file[i-1]; }
#endif /* DOXYGEN_SHOULD_SKIP_THIS */

/// \brief Return <tt>i</tt>-th parameter read using keyword \b Auxfile.
/// \details Here we have at most 10 integer extra parameters that can be used for any auxiliary file names.
/// Default value for <tt>i</tt> is <tt>1</tt>
    string getAuxFile(size_t i=1) const { return _aux_file[i-1]; }

/// \brief Return expression (to be parsed, function of <tt>x</tt>, <tt>y</tt>, <tt>z</tt>,
/// <tt>t</tt>) for density function.
    string getDensity() const { return _mp._density; }

/// \brief Return expression (to be parsed, function of <tt>x</tt>, <tt>y</tt>, <tt>z</tt>,
/// <tt>t</tt>) for electric conductivity.
    string getElectricConductivity() const { return _mp._electric_cond; }

/// \brief Return expression (to be parsed, function of <tt>x</tt>, <tt>y</tt>, <tt>z</tt>,
/// <tt>t</tt>) for electric permittivity.
    string getElectricPermittivity() const { return _mp._electric_perm; }

/// \brief Return expression (to be parsed, function of <tt>x</tt>, <tt>y</tt>, <tt>z</tt>,
/// <tt>t</tt>) for magnetic permeability.
    string getMagneticPermeability() const { return _mp._magnetic_perm; }

/// \brief Return expression (to be parsed, function of <tt>x</tt>, <tt>y</tt>, <tt>z</tt>,
/// <tt>t</tt>) for Poisson ratio.
    string getPoissonRatio() const { return _mp._poisson; }

/// \brief Return expression (to be parsed, function of <tt>x</tt>, <tt>y</tt>, <tt>z</tt>,
/// <tt>t</tt>) for thermal conductivity.
    string getThermalConductivity() const { return _mp._thermal_cond; }

/// \brief Return expression (to be parsed, function of <tt>x</tt>, <tt>y</tt>, <tt>z</tt>,
/// <tt>t</tt>) for density * specific heat.
    string getRhoCp() const { return _mp._rho_cp; }

/// \brief Return expression (to be parsed, function of <tt>x</tt>, <tt>y</tt>, <tt>z</tt>,
/// <tt>t</tt>) for viscosity.
    string getViscosity() const { return _mp._visc; }

/// \brief Return expression (to be parsed, function of <tt>x</tt>, <tt>y</tt>, <tt>z</tt>,
/// <tt>t</tt>) for Young's modulus.
    string getYoungModulus() const { return _mp._young; }

    friend class XMLParser;

private:

   int _bc, _sf, _bf, _ini, _data;
   size_t _nb_par, _nb_steps, _nb_iter, _nb_aux_files, _nb_mesh_files;
   int _verbose, _save, _plot, _output;
   size_t _nb_plot_files, _nb_complex_par, _nb_string_par, _nb_data_files;
   real_t _tolerance, _time_step, _max_time;
   string _file, _project, _domain_file, _init_file, _restart_file;
   string _bc_file, _bf_file, _sf_file, _save_file;
   size_t _nb_mat, _nb_int_par, _nb_real_par, _nb_point_double_par;
   size_t _nmat[MAX_NB_MATERIALS];
   int _int_par[MAX_NB_PAR];
   real_t _real_array[MAX_NB_PAR][MAX_ARRAY_SIZE];
   real_t _real_par[MAX_NB_PAR];
   Point<real_t> _point_double_par[MAX_NB_PAR];
   complex<real_t> _complex_par[MAX_NB_PAR];
   string _mat[MAX_NB_MATERIALS], _mesh_file[MAX_NB_PAR], _aux_file[MAX_NB_PAR];
   string _data_file[MAX_NB_PAR], _plot_file[MAX_NB_PAR], _string_par[MAX_NB_PAR];
   vector<string> _param_label, _param_value, _param_ext;
   vector<string> _array_label, _array_ext;
   vector<size_t> _array_size;
   vector<real_t> _array_value[MAX_ARRAY_SIZE];
   mat_prop _mp;
   void init();
};

/*! @} End of Doxygen Groups */
} /* namespace OFELI */

#endif
