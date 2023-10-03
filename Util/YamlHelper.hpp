/**
 * @file YamlHelper.hpp
 * @author Wen Yan (wenyan4work@gmail.com)
 * @brief  type deduction and error message
 * @version 0.1
 * @date 2020-01-03
 *
 * @copyright Copyright (c) 2020
 *
 */
#ifndef YAMLHELPER_HPP_
#define YAMLHELPER_HPP_

#include "Logger.hpp"

#include <iostream>
#include <yaml-cpp/yaml.h>

/**
 * @brief convert variable name to string
 *
 */
#define VARNAME(a) #a

/**
 * @brief Read one variable from
 *
 * @tparam T  type of variable
 * @param config yaml config file
 * @param name name to look in yaml
 * @param variable [out] value found in yaml
 * @param err extra error message
 * @param optional if this variable is optional
 */
template <typename T>
void readConfig(const YAML::Node &config, const std::string &name, T &variable, const std::string &err,
                const bool optional = false) {
    if (config[name]) {
        variable = config[name].as<T>();
    } else {
        if (!optional) {
            spdlog::critical("Required parameter " + name + " in input yaml file not found");
            std::exit(1);
        }
        spdlog::warn("Optional " + name + " in input yaml file not found, using default");
    }
}

/**
 * @brief fixed-size array version of readConfig
 *
 * @tparam T
 * @param config
 * @param name
 * @param variable
 * @param err
 * @param optional
 */
template <typename T>
void readConfig(const YAML::Node &config, const std::string &name, T variable[], int dim, const std::string &err,
                const bool optional = false) {
    if (config[name]) {
        YAML::Node seq = config[name];
        if (seq.size() != dim) {
            spdlog::critical("Expecting {} elements in {} in input yaml file. {}. ", dim, name, err);
            std::exit(1);
        }
        for (int i = 0; i < dim; i++) {
            variable[i] = seq[i].as<T>();
        }
    } else {
        if (!optional) {
            spdlog::critical("Required parameter " + name + " in input yaml file not found");
            std::exit(1);
        }
        spdlog::warn("Optional " + name + " in input yaml file not found, using default");
    }
}

/**
 * @brief fixed-size array version of readConfig
 *
 * @tparam T
 * @param config
 * @param name
 * @param variable
 * @param err
 * @param optional
 */
template <typename T>
void readVariableConfig(const YAML::Node &config, const std::string &name, T variable[], int dim, const std::string &err,
                const bool optional = false) {
    if (config[name]) {
        YAML::Node seq = config[name];
        if (seq.IsSequence()) {
            if (seq.size() != dim) {
                spdlog::critical("Expecting {} elements in {} in input yaml file. {}. ", dim, name, err);
                std::exit(1);
            }
            for (int i = 0; i < dim; i++) {
                variable[i] = seq[i].as<T>();
            }
        } else if (seq.IsScalar()){
            variable[0] = seq.as<T>();
            for (int i = 1; i < dim; i++) {
                variable[i] = variable[0];
            }
        } else{
                spdlog::critical("Node type was not what was expected for {} in input yaml file. {}. ", dim, name, err);
                std::exit(1);
            }
    } else {
        if (!optional) {
            spdlog::critical("Required parameter " + name + " in input yaml file not found");
            std::exit(1);
        }
        spdlog::warn("Optional " + name + " in input yaml file not found, using default");
    }
}

/**
 * @brief std::vector version of
 *
 * @tparam T
 * @param config
 * @param name
 * @param variable
 * @param err
 * @param optional
 */
template <typename T>
void readConfig(const YAML::Node &config, const std::string &name, std::vector<T> &variable, const std::string &err,
                const bool optional = false) {
    if (config[name]) {
        YAML::Node seq = config[name];
        const int dim = seq.size();
        variable.resize(dim);
        for (int i = 0; i < dim; i++) {
            variable[i] = seq[i].as<T>();
        }
    } else {
        if (!optional) {
            spdlog::critical("Required parameter " + name + " in input yaml file not found");
            std::exit(1);
        }
        spdlog::warn("Optional " + name + " in input yaml file not found, using default");
    }
}

#endif