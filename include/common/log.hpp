#ifndef NONLINEAR_CONTROLS_COMMON_LOG_HPP
#define NONLINEAR_CONTROLS_COMMON_LOG_HPP

#include "common/utils.hpp"
#include <iostream>

namespace nonlinear_controls {

static const struct TerminalColors {
  std::string GREEN = "\033[01;32m";
  std::string NC = "\033[0m"; // No Color
  std::string BLACK = "\033[01;30m";
  std::string RED = "\033[01;31m";
  std::string YELLOW = "\033[01;33m";
  std::string BLUE = "\033[01;34m";
  std::string MAGENTA = "\033[01;35m";
  std::string CYAN = "\033[01;36m";
  std::string WHITE = "\033[0;37m";
  std::string Reset = "\033[0m";
} tColors;

template<typename T> std::string Cat(T arg) {
/* http://www.cplusplus.com/forum/beginner/235373/ */
  std::stringstream ss;
  ss << arg;
  return ss.str();
}

template<typename T, typename... Args>
std::string Cat(T current, Args... args) {
  std::string result;
  result += Cat(current);
  result += Cat((args)...);
  return result;
}

class Logger {

public:
  Logger() = default;

  ~Logger() = default;

  static std::string get_curr_time_as_string() {
    return std::to_string(static_cast<double>(utils::get_current_time()) * 1.0e-6) + "s";
  }

  static void WARN(const std::string &msg) {
    std::string s;
    s += std::string("[WARN] [") + get_curr_time_as_string() +
        std::string("]: ") + msg;
    std::cout << tColors.YELLOW << s << tColors.NC << std::endl;
  }

  static void ERROR(const std::string &msg) {
    std::string s;
    s += std::string("[ERROR] [") + get_curr_time_as_string() +
        std::string("]: ") + msg;
    std::cout << tColors.RED + s + tColors.NC << std::endl;
  }
  static void SUCCESS(const std::string &msg) {
    std::string s;
    s += std::string("[SUCCESS] [") +
        get_curr_time_as_string() + std::string("]: ") + msg;
    std::cout << tColors.GREEN + s + tColors.NC << std::endl;
  }

  static void STATUS(const std::string &msg) {
    std::string s;
    s += std::string("[STATUS] [") + get_curr_time_as_string() +
        std::string("]: ") + msg;
    std::cout << tColors.CYAN + s + tColors.NC << std::endl;
  }

  static void INFO(const std::string &msg) {
    std::string s;
    s += std::string("[INFO] [") + get_curr_time_as_string() +
        std::string("]: ") + msg;
    std::cout << tColors.WHITE + s + tColors.NC << std::endl;
  }
};

}

#endif //NONLINEAR_CONTROLS_COMMON_LOG_HPP
