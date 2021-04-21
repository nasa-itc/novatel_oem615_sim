/* Copyright (C) 2015 - 2019 National Aeronautics and Space Administration. All Foreign Rights are Reserved to the U.S. Government.

   This software is provided "as is" without any warranty of any, kind either express, implied, or statutory, including, but not
   limited to, any warranty that the software will conform to, specifications any implied warranties of merchantability, fitness
   for a particular purpose, and freedom from infringement, and any warranty that the documentation will conform to the program, or
   any warranty that the software will be error free.

   In no event shall NASA be liable for any damages, including, but not limited to direct, indirect, special or consequential damages,
   arising out of, resulting from, or in any way connected with the software or its documentation.  Whether or not based upon warranty,
   contract, tort or otherwise, and whether or not loss was sustained from, or arose out of the results of, or use of, the software,
   documentation or services provided hereunder

   ITC Team
   NASA IV&V
   ivv-itc@lists.nasa.gov
*/

#include <iomanip>
#include <string>
#include <cmath>

#include <ItcLogger/Logger.hpp>

#include <sim_coordinate_transformations.hpp>
#include <gps_sim_data_point.hpp>

namespace Nos3
{
    extern ItcLogger::Logger *sim_logger;

    /*************************************************************************
     * Constructors
     *************************************************************************/

    GPSSimDataPoint::GPSSimDataPoint(double abs_time, int16_t gps_week, int32_t gps_sec_week, double gps_frac_sec, 
        const std::vector<double>& ECEF, const std::vector<double>& ECEF_vel, const std::vector<double>& ECI, const std::vector<double>& ECI_vel) : 
        _not_parsed(false), _abs_time(abs_time), _gps_week(gps_week), _gps_sec_week(gps_sec_week), _gps_frac_sec(gps_frac_sec), _ECEF(ECEF), _ECEF_vel(ECEF_vel), _ECI(ECI), _ECI_vel(ECI_vel)
    {
    }

    GPSSimDataPoint::GPSSimDataPoint(int16_t spacecraft, int16_t gps, const boost::shared_ptr<Sim42DataPoint> dp) : 
        _sc(spacecraft), _gps(gps), _dp(*dp), _not_parsed(true) 
    {
        _ECEF.resize(3); _ECEF_vel.resize(3); _ECI.resize(3); _ECI_vel.resize(3);
        sim_logger->trace("GPSSimDataPoint::GPSSimDataPoint:  Created instance using _sc=%d, _gps=%d, _dp=%s", 
            _sc, _gps, _dp.to_string().c_str());
    }

    /*************************************************************************
     * Mutators
     *************************************************************************/

    void GPSSimDataPoint::do_parsing(void) const
    {
        std::ostringstream MatchString;
        MatchString << "SC[" << _sc << "].AC.GPS[" << _gps << "].";
        size_t MSsize = MatchString.str().size();

        _not_parsed = false;
        
        std::vector<std::string> lines = _dp.get_lines();
        
        try {
            for (int i = 0; i < lines.size(); i++) {
                if (lines[i].compare(0, MSsize, MatchString.str()) == 0) { // e.g. SC[0].AC.GPS[0]
                    sim_logger->trace("GPSSimDataPoint::do_parsing:  Found a string with the correct prefix = %s.  String:  %s", MatchString.str().c_str(), lines[i].c_str());
                    // GPS. Rollover, Week, Sec, PosN, VelN, PosW, VelW, Lng, Lat, Alt
                    if (lines[i].compare(MSsize, 9, "Rollover ") == 0) {
                        size_t found = lines[i].find_first_of("=");
                        _gps_rollover = std::stoi(lines[i].substr(found+1, lines[i].size() - found - 1));
                        sim_logger->trace("GPSSimDataPoint::do_parsing:  Parsed Rollover.  = found at %d, rhs=%s, _gps_rollover=%d", 
                            found, lines[i].substr(found+1, lines[i].size() - found - 1).c_str(), _gps_rollover);
                    } else
                    if (lines[i].compare(MSsize, 5, "Week ") == 0) {
                        size_t found = lines[i].find_first_of("=");
                        _gps_week = std::stoi(lines[i].substr(found+1, lines[i].size() - found - 1));
                        sim_logger->trace("GPSSimDataPoint::do_parsing:  Parsed Week.  = found at %d, rhs=%s, _gps_week=%d", 
                            found, lines[i].substr(found+1, lines[i].size() - found - 1).c_str(), _gps_week);
                    } else
                    if (lines[i].compare(MSsize, 4, "Sec ") == 0) {
                        size_t found = lines[i].find_first_of("=");
                        double seconds = std::stod(lines[i].substr(found+1, lines[i].size() - found - 1));
                        _gps_sec_week = (int32_t)seconds;
                        _gps_frac_sec = seconds - (double)_gps_sec_week;
                        sim_logger->trace("GPSSimDataPoint::do_parsing:  Parsed Sec.  = found at %d, rhs=%s, seconds=%f, _gps_sec_week=%d, _gps_frac_sec=%f", 
                            found, lines[i].substr(found+1, lines[i].size() - found - 1).c_str(), seconds, _gps_sec_week, _gps_frac_sec);
                    } else
                    if (lines[i].compare(MSsize, 5, "PosN ") == 0) {
                        size_t found = lines[i].find_first_of("=");
                        parse_double_vector(lines[i].substr(found+1, lines[i].size() - found - 1), _ECI);
                        sim_logger->trace("GPSSimDataPoint::do_parsing:  Parsed PosN.  = found at %d, rhs=%s, _ECI=%f/%f/%f", 
                            found, lines[i].substr(found+1, lines[i].size() - found - 1).c_str(), _ECI[0], _ECI[1], _ECI[2]);
                    } else
                    if (lines[i].compare(MSsize, 5, "VelN ") == 0) {
                        size_t found = lines[i].find_first_of("=");
                        parse_double_vector(lines[i].substr(found+1, lines[i].size() - found - 1), _ECI_vel);
                        sim_logger->trace("GPSSimDataPoint::do_parsing:  Parsed VelN.  = found at %d, rhs=%s, _ECI_vel=%f/%f/%f", 
                            found, lines[i].substr(found+1, lines[i].size() - found - 1).c_str(), _ECI_vel[0], _ECI_vel[1], _ECI_vel[2]);
                    } else
                    if (lines[i].compare(MSsize, 5, "PosW ") == 0) {
                        size_t found = lines[i].find_first_of("=");
                        parse_double_vector(lines[i].substr(found+1, lines[i].size() - found - 1), _ECEF);
                        sim_logger->trace("GPSSimDataPoint::do_parsing:  Parsed PosW.  = found at %d, rhs=%s, _ECEF=%f/%f/%f", 
                            found, lines[i].substr(found+1, lines[i].size() - found - 1).c_str(), _ECEF[0], _ECEF[1], _ECEF[2]);
                    } else
                    if (lines[i].compare(MSsize, 5, "VelW ") == 0) {
                        size_t found = lines[i].find_first_of("=");
                        parse_double_vector(lines[i].substr(found+1, lines[i].size() - found - 1), _ECEF_vel);
                        sim_logger->trace("GPSSimDataPoint::do_parsing:  Parsed VelW.  = found at %d, rhs=%s, _ECEF_vel=%f/%f/%f", 
                            found, lines[i].substr(found+1, lines[i].size() - found - 1).c_str(), _ECEF_vel[0], _ECEF_vel[1], _ECEF_vel[2]);
                    } else
                    if (lines[i].compare(MSsize, 4, "Lng ") == 0) {
                        size_t found = lines[i].find_first_of("=");
                        _gps_lng = std::stod(lines[i].substr(found+1, lines[i].size() - found - 1)) * 180.0 / M_PI;
                        sim_logger->trace("GPSSimDataPoint::do_parsing:  Parsed Lng.  = found at %d, rhs=%s, _gps_lng=%f", 
                            found, lines[i].substr(found+1, lines[i].size() - found - 1).c_str(), _gps_lng);
                    } else
                    if (lines[i].compare(MSsize, 4, "Lat ") == 0) {
                        size_t found = lines[i].find_first_of("=");
                        _gps_lat = std::stod(lines[i].substr(found+1, lines[i].size() - found - 1)) * 180.0 / M_PI;
                        sim_logger->trace("GPSSimDataPoint::do_parsing:  Parsed Lat.  = found at %d, rhs=%s, _gps_lat=%f", 
                            found, lines[i].substr(found+1, lines[i].size() - found - 1).c_str(), _gps_lat);
                    } else
                    if (lines[i].compare(MSsize, 4, "Alt ") == 0) {
                        size_t found = lines[i].find_first_of("=");
                        _gps_alt = std::stod(lines[i].substr(found+1, lines[i].size() - found - 1));
                        sim_logger->trace("GPSSimDataPoint::do_parsing:  Parsed Alt.  = found at %d, rhs=%s, _gps_alt=%f", 
                            found, lines[i].substr(found+1, lines[i].size() - found - 1).c_str(), _gps_alt);
                    }
                }
            }
        } catch(const std::exception& e) {
            sim_logger->error("GPSSimDataPoint::do_parsing:  Parsing exception:  %s", e.what());
        }
        
        double JD;
        SimCoordinateTransformations::GpsTimeToJD(_gps_rollover, _gps_week, ((double)_gps_sec_week) + _gps_frac_sec, JD);
        _abs_time = SimCoordinateTransformations::JDToAbsTime(JD);
        
        sim_logger->debug("GPSSimDataPoint::do_parsing:  Parsed data point:\n%s", to_string().c_str());

    }

    /*************************************************************************
     * Accessors
     *************************************************************************/

    std::string GPSSimDataPoint::to_formatted_string(void) const
    {
        parse_data_point();
        
        std::stringstream ss;

        ss << std::fixed << std::setprecision(4) << std::setfill(' ');
        ss << "GPS Data Point: " << std::endl;
        ss << "  Absolute Time                    : " << std::setw(15) << _abs_time << std::endl;
        ss << "  GPS Rollover, Week, Second, Fractional Second: "
           << std::setw(6) << _gps_rollover << ","
           << std::setw(6) << _gps_week << ","
           << std::setw(7) << _gps_sec_week << ","
           << std::setw(7) << _gps_frac_sec << std::endl;
        ss << std::setprecision(2);
        ss << "  ECEF                                         : "
           << std::setw(12) << _ECEF[0] << ","
           << std::setw(12) << _ECEF[1] << ","
           << std::setw(12) << _ECEF[2] << std::endl;
        ss << "  ECEF Velocity                                : "
           << std::setw(12) << _ECEF_vel[0] << ","
           << std::setw(12) << _ECEF_vel[1] << ","
           << std::setw(12) << _ECEF_vel[2] << std::endl;
        ss << "  ECI                                          : "
           << std::setw(12) << _ECI[0] << ","
           << std::setw(12) << _ECI[1] << ","
           << std::setw(12) << _ECI[2] << std::endl;
        ss << "  ECI Velocity                                 : "
           << std::setw(12) << _ECI_vel[0] << ","
           << std::setw(12) << _ECI_vel[1] << ","
           << std::setw(12) << _ECI_vel[2] << std::endl;
        ss << "  Geodetic Lat/Lng/Alt(m above WGS-84)         : "
           << std::setw(12) << _gps_lat << ","
           << std::setw(12) << _gps_lng << ","
           << std::setw(12) << _gps_alt << std::endl;

        return ss.str();
    }

    std::string GPSSimDataPoint::to_string(void) const
    {
        parse_data_point();
        
        std::stringstream ss;

        ss << std::fixed << std::setfill(' ');
        ss << "GPS Data Point: ";
        ss << std::setprecision(std::numeric_limits<double>::digits10); // Full double precision
        ss << " AbsTime: " << _abs_time;
        ss << std::setprecision(std::numeric_limits<int32_t>::digits10); // Full int32_t precision
        ss << " GPS Time: "
           << _gps_rollover << "/"
           << _gps_week << "/"
           << _gps_sec_week << "/";
        ss << std::setprecision(std::numeric_limits<double>::digits10); // Full double precision
        ss << _gps_frac_sec ;
        ss << " ECEF: "
           << _ECEF[0] << ","
           << _ECEF[1] << ","
           << _ECEF[2] ;
        ss << " ECEF Velocity"
           << _ECEF_vel[0] << ","
           << _ECEF_vel[1] << ","
           << _ECEF_vel[2] ;
        ss << " ECI "
           << _ECI[0] << ","
           << _ECI[1] << ","
           << _ECI[2] ;
        ss << " ECI Velocity"
           << _ECI_vel[0] << ","
           << _ECI_vel[1] << ","
           << _ECI_vel[2] ;
        ss << "  Geodetic Lat/Lng/Alt(m above WGS-84): "
           << _gps_lat << ","
           << _gps_lng << ","
           << _gps_alt ;

        return ss.str();
    }
    
}
