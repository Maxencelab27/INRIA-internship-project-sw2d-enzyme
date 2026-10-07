#pragma once

#include <string>

class swIO;

class swPostProcessParam
{
public :
    swPostProcessParam(swIO *io);
    ~swPostProcessParam() = default;

    //-- General methods
    bool checkConsistency();
    bool readInputFile();
    bool writeLog();
    //== General methods

    //-- Getters
    bool getComputeMax();
    bool getFirstSubmersion();
    double getHmin();
    std::string getInputFormat();
    std::string getInputMax();
    std::string getInputSubmersion();
    std::string getSubmersionDurations();
    //== Getters

private:
    //-- Setters
    bool setcomputeMax(const std::string& text);
    bool setFirstSubmersion(const std::string& text);
    bool setHmin(const std::string& text);
    bool setInputFormat(const std::string& text);
    bool setInputMax(const std::string& text);
    bool setInputSubmersion(const std::string& text);
    bool setSubmersionDurations(const std::string& text);
    //== Setters

    //-- Attributes
    swIO *m_io = nullptr;

    bool compute_max = false;
    bool first_submersion = false;
    double h_min = 0.01;
    std::string input_format = "txt";
    std::string input_max = "all";
    std::string input_submersion = "all";
    std::string submersion_durations = "none";
    //== Attributes

};
