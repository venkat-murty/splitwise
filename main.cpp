
#include "csv.hpp"
#include <algorithm>
#include <iostream>
#include <numeric>
#include <regex>
#include <set>
#include <stdio.h>
#include <string>

struct Row : public std::vector<std::string> {
  int line;
  Row(int idx) : std::vector<std::string>{}, line(idx) {}
};
using Data = std::vector<Row>;

Row nil{0};
Row &curr = nil;

std::vector<std::string> dates;
std::map<int, std::string> peeps;

void check(int line, const char *msg, bool value) {
  if (value)
    return;

  std::cerr << "Code Line " << line << " " << msg << std::endl;
  std::cerr << "CSV Line " << curr.line << std::endl;
  std::for_each(curr.begin(), curr.end(),
                [](const auto &l) { std::cerr << l << "\t"; });
  std::cerr << std::endl;

  abort();
}

#define CHECK(expr, msg) check(__LINE__, msg, expr)

std::string ltrim(const std::string &s) {
  return std::regex_replace(s, std::regex("^\\s+"), std::string(""));
}

std::string rtrim(const std::string &s) {
  return std::regex_replace(s, std::regex("\\s+$"), std::string(""));
}

std::string trim(const std::string &s) { return ltrim(rtrim(s)); }

bool isascii(const std::string &s) {
  for (char c : s) {
    if (isspace(c) || isalpha(c) || isdigit(c))
      continue;
    return false;
  }
  return true;
}

bool isdate(std::string const &d) {
  for (const auto &x : dates)
    if (x == d)
      return true;
  return false;
}

bool ispeep(std::string const &p) {
  for (const auto &x : peeps)
    if (x.second == p)
      return true;
  return false;
}

bool ispeep(int idx) { return peeps.find(idx) != peeps.end(); }

Data p_data() {
  std::ifstream file{"expense.csv"};
  csv::Reader csv{file}; // open a std::istream

  Data result;

  int idx = 1;
  for (auto &&row : csv) // row is a csv::Row object
  {
    Row r{idx++};
    for (auto &&column : row) // column is a std::string. To get another type,
                              // use Row::Range. For example:
                              // for(auto && column: row.range<int>())
    {
      r.push_back(trim(column));
    }
    result.emplace_back(std::move(r));
  }
  return result;
}

void p_dates(Data const &data) {
  std::vector<std::string> &result = dates;
  for (auto &row : data) {
    curr = row;

    if (row.at(0) == "Dates") {
      for (auto &col : row) {
        if (!(col.empty() || col == "Dates")) {
          CHECK(isdate(col) == false, "Duplicate date");
          result.push_back(col);
        }
      }
      break;
    }
  }
  curr = nil;

  //
  //
  std::cout << "Dates ";
  std::for_each(result.begin(), result.end(),
                [](const auto &d) { std::cout << "'" << d << "' "; });
  std::cout << std::endl;
}

void p_peeps(Data const &data) {
  std::map<int, std::string> &result = peeps;

  for (auto &row : data) {
    curr = row;

    if (row.at(0) == "Peeps") {
      for (int i = 1; i < row.size(); ++i) {
        auto &col = row[i];
        if (!col.empty()) {
          CHECK(ispeep(col) == false, "Duplicate Peeps");
          result[i] = col;
        }
      }
      break;
    }
  }
  curr = nil;

  std::cout << "Peeps ";
  std::for_each(result.begin(), result.end(),
                [](const auto &d) { std::cout << "'" << d.second << "' "; });
  std::cout << std::endl;
}

void p_echo(Data const &data) {
  auto d = "";
  for (auto &row : data) // row is a csv::Row object
  {
    for (auto &col : row) {
      std::cout << d << col << d << "\t";
    }
    std::cout << std::endl;
  }
}

bool isnumber(std::string const &n) {
  for (char c : n) {
    if (c == '.' || isdigit(c))
      continue;
    else {
      std::cerr << "Expecting a number and got " << n << std::endl;
      return false;
    }
  }
  return true;
}

Data subset(Data const &data, std::string token, std::string next_token) {
  Data result;

  int r = 0;
  for (; r < data.size(); ++r) {
    if (data.at(r).empty())
      continue;
    if (data.at(r).at(0) == token)
      break;
  }
  CHECK(data.size() > r, "Token not found in the input");
  CHECK(data.at(r).at(0) == token, "Token not found");

  if (data.at(r).at(1) == "Amount")
    ++r;

  for (; r < data.size(); ++r) {
    if (data.at(r).empty())
      continue;
    if (data.at(r).at(0) == next_token || data.at(r).at(0) == "Notes")
      return result;

    if (!data.at(r).at(1).empty()) {
      result.push_back(data.at(r));
    }
  }
  return result;
}

// bookie => {person => payment}
std::map<std::string, std::map<std::string, double>>
p_bookies(Data const &data) {
  std::map<std::string, std::map<std::string, double>> result;

  for (auto &row : data) {
    curr = row;

    auto bookie = row.at(1);
    CHECK(ispeep(bookie), "Bookie is not valid");
    CHECK(result.find(bookie) == result.end(), "Multiple line per bookie");

    result[bookie] = {};

    for (int i = 2; i < row.size(); ++i) {
      auto &col = row[i];
      if (!col.empty()) {
        CHECK(isnumber(col), "Amount should be a number");
        CHECK(ispeep(i), "Invalid data");

        result[bookie][peeps.at(i)] = std::stod(col);
      }
    }
  }
  curr = nil;

  std::cout << "Bookies and initialize deposit " << std::endl;
  for (auto &b : result) {
    std::cout << "\t" << b.first << std::endl;
    for (auto &p : b.second) {
      std::cout << "\t\t" << p.first << "\t= " << p.second << std::endl;
    }
  }
  return result;
}

// Date => {peep => wt}
std::map<std::string, std::map<std::string, double>>
p_attendence(Data const &data) {
  std::map<std::string, std::map<std::string, double>> result;

  for (auto &row : data) {
    curr = row;

    std::string dt = row.at(1);
    CHECK(isdate(dt), "Ivalid date in attendence");

    for (int i = 2; i < row.size(); ++i) {
      auto &col = row.at(i);
      if (!col.empty()) {
        CHECK(isnumber(col), "Amount should be a number");
        CHECK(ispeep(i), "Invalid assignment");

        result[dt][peeps.at(i)] = std::stod(col);
      }
    }
  }
  curr = nil;

  std::cout << "Number of Peeps" << std::endl;
  std::for_each(result.begin(), result.end(), [](auto const &date_kv) {
    std::cout << "\t" << date_kv.first << std::endl;
    std::for_each(
        date_kv.second.begin(), date_kv.second.end(), [](auto const &kv) {
          std::cout << "\t\t" << kv.first << "\t= " << kv.second << std::endl;
        });
  });

  return result;
}

using Expense =
    std::tuple<std::string, double, std::map<std::string, double>, std::string>;

// {person, amount, {person => wt}, desc}
std::vector<Expense> p_assigned_expense(Data const &data) {
  std::vector<Expense> result;

  for (auto &row : data) {
    curr = row;
    CHECK(isnumber(row.at(1)), "Invalid amount in assigned expenses");
    CHECK(ispeep(row.at(0)), "Invalid person in assigned expenses");

    auto payer = row.at(0);
    auto amount = std::stod(row.at(1));
    auto desc = row.at(2);

    Expense split{payer, amount, {}, desc};
    for (int i = 3; i < row.size(); ++i) {
      auto &col = row.at(i);
      if (!col.empty()) {
        CHECK(isnumber(col), "Invalid wt in assigned expenses");
        CHECK(ispeep(i), "Invalid assignment in assigned expenses");

        std::get<2>(split)[peeps.at(i)] = std::stod(col);
      }
    }
    if (!std::get<2>(split).empty()) {
      result.emplace_back(std::move(split));
    }
  }
  curr = nil;

  //
  std::cout << "Assigned Amount" << std::endl;
  std::for_each(result.begin(), result.end(), [](auto const &elem) {
    std::cout << "\t" << std::get<0>(elem) << "\t" << std::get<1>(elem) << "\t"
              << std::get<3>(elem) << std::endl;

    std::for_each(
        std::get<2>(elem).begin(), std::get<2>(elem).end(), [](auto const &kv) {
          std::cout << "\t\t" << kv.first << "\t= " << kv.second << std::endl;
        });
  });
  return result;
}

int date_idx(std::string const &d) {
  for (int i = 0; i < dates.size(); ++i)
    if (dates.at(i) == d)
      return i;
  return -1;
}

std::vector<Expense> p_shared_expenses(
    Data const &data,
    std::map<std::string, std::map<std::string, double>> const &attendence) {

  std::vector<Expense> result;
  for (auto &row : data) {
    curr = row;

    CHECK(row.size() >= 3, "Invalid shared expenses");

    CHECK(ispeep(row.at(0)), "Expecting a person in first column");
    CHECK(isnumber(row.at(1)), "Expecting an amount in second column");
    CHECK(isdate(row.at(2)), "Expecting a date in third column");

    auto &payer = row.at(0);
    auto amount = std::stod(row.at(1));
    auto &start = row.at(2);
    auto &end = row.size() >= 4 && !row.at(3).empty() ? row.at(3) : row.at(2);
    std::string desc = row.size() >= 5 ? row.at(4) : "";

    CHECK(isdate(end), "Expecting a date in fourth column");

    int s = date_idx(start);
    int e = date_idx(end);

    CHECK(s != -1, "Program Bug");
    CHECK(e != -1, "Program Bug");
    CHECK(s <= e, "Invalid date range");

    double days = 1.0 + (e - s);

    for (int d = s; d <= e; ++d) {
      std::string const &date = dates.at(d);
      CHECK(attendence.find(date) != attendence.end(),
            "Attendence matrix not complete");

      Expense split{payer, amount / days, attendence.find(date)->second, desc};
      result.emplace_back(std::move(split));
    }
  }

  curr = nil;

  //
  std::cout << "Shared Amount" << std::endl;
  std::for_each(result.begin(), result.end(), [](auto const &elem) {
    std::cout << "\t" << std::get<0>(elem) << "\t" << std::get<1>(elem) << "\t"
              << std::get<3>(elem) << std::endl;

    std::for_each(
        std::get<2>(elem).begin(), std::get<2>(elem).end(), [](auto const &kv) {
          std::cout << "\t\t" << kv.first << "\t= " << kv.second << std::endl;
        });
  });
  return result;
}

std::vector<Expense> p_initial_expenses(
    std::map<std::string, std::map<std::string, double>> const &bookies) {
  std::vector<Expense> result;
  for (auto &elem : bookies) {
    for (auto &e : elem.second) {
      Expense split{e.first, e.second, {{elem.first, 1.0}}, "Initial Deposit"};
      result.emplace_back(std::move(split));
    }
  }

  std::cout << "Initial Amount" << std::endl;
  std::for_each(result.begin(), result.end(), [](auto const &elem) {
    std::cout << "\t" << std::get<0>(elem) << "\t" << std::get<1>(elem) << "\t"
              << std::get<3>(elem) << std::endl;

    std::for_each(
        std::get<2>(elem).begin(), std::get<2>(elem).end(), [](auto const &kv) {
          std::cout << "\t\t" << kv.first << "\t= " << kv.second << std::endl;
        });
  });
  return result;
}

void adjust(std::map<std::string, double> &amounts,
            std::vector<Expense> const &expenses) {
  for (Expense const &expense : expenses) {
    auto const &wts = std::get<std::map<std::string, double>>(expense);

    auto total = std::accumulate(
        wts.begin(), wts.end(), 0.0,
        [](double sum, auto const &wt) { return sum + wt.second; });

    auto payment = std::get<double>(expense);

    amounts.find(std::get<0>(expense))->second += payment;

    std::for_each(
        wts.begin(), wts.end(), [&amounts, payment, total](auto const &wt) {
          amounts.find(wt.first)->second -= payment * wt.second / total;
        });
  }
}

void echo(std::map<std::string, double> &amounts) {
  double t = 0.0;
  std::cout << "Amounts" << std::endl;
  std::for_each(amounts.begin(), amounts.end(), [&t](auto const &e) {
    std::cout << "\t" << e.first << "\t= " << e.second << std::endl;
    t += e.second;
  });
  std::cout << "\tSum\t= " << t << std::endl;
}

// Bookies and expenses.
std::tuple<std::vector<std::string>, std::map<std::string, double>> expenses() {
  auto data = p_data();
  p_echo(data);
  std::cout << std::endl << std::endl;

  p_dates(data);
  p_peeps(data);

  std::map<std::string, std::map<std::string, double>> bookies =
      p_bookies(subset(data, "Bookies", "Number of Peeps"));

  // Date => {peep => wt}
  std::map<std::string, std::map<std::string, double>> attendence =
      p_attendence(subset(data, "Number of Peeps", "Assigned Expenses"));

  // {person, amount, {person => wt}, desc}
  std::vector<Expense> assigned_expenses =
      p_assigned_expense(subset(data, "Assigned Expenses", "Shared Expenses"));

  std::vector<Expense> shared_expenses =
      p_shared_expenses(subset(data, "Shared Expenses", "Notes"), attendence);

  std::vector<Expense> initial_expenses = p_initial_expenses(bookies);

  auto amounts = std::accumulate(
      peeps.begin(), peeps.end(), std::map<std::string, double>{},
      [](std::map<std::string, double> result, auto const &p) {
        result[p.second] = 0.0;
        return result;
      });

  adjust(amounts, shared_expenses);
  adjust(amounts, assigned_expenses);
  adjust(amounts, initial_expenses);

  echo(amounts);

  return {std::accumulate(
              bookies.begin(), bookies.end(), std::vector<std::string>{},
              [](std::vector<std::string> result, auto const &elem) {
                result.push_back(elem.first);
                return result;
              }),
          amounts};
}

#include "lp_lib.h"

int main() {
  const auto [bookies, amounts] = expenses();

  std::vector<std::pair<std::string, std::string>> payments;
  for (auto const &b : bookies) {
    for (auto const &a : amounts) {
      if (a.first != b) {
        payments.push_back({a.first, b});
        payments.push_back({b, a.first});
      }
    }
  }

  lprec *lp;

  int *colno = NULL;
  REAL *row = NULL;

  const int Ncol = payments.size() * 2;

  lp = make_lp(0, Ncol);

  std::map<std::string, std::vector<int>> inflow;
  std::map<std::string, std::vector<int>> outflow;

  std::set<int> binary;

  // Payment variables.
  for (int i = 0; i < payments.size(); ++i) {
    std::string a = payments.at(i).first + "->" + payments.at(i).second;
    std::string b = a + "_b";

    set_col_name(lp, i * 2 + 1, const_cast<char *>(a.c_str()));
    set_bounds(lp, i * 2 + 1, 0.0, get_infinite(lp));

    outflow[payments.at(i).first].push_back(i * 2 + 1);
    inflow[payments.at(i).second].push_back(i * 2 + 1);

    set_col_name(lp, i * 2 + 2, const_cast<char *>(b.c_str()));
    set_binary(lp, i * 2 + 2, TRUE);

    binary.insert(i * 2 + 2);
  }
  // allocate
  colno = (int *)malloc(Ncol * sizeof(int));
  row = (REAL *)malloc(Ncol * sizeof(REAL));

  set_add_rowmode(lp, TRUE);

  // Binary variable constraint
  for (int i = 0; i < payments.size(); ++i) {
    int vcol = i * 2 + 1;
    int bcol = i * 2 + 2;

    int j = 0;

    const REAL SMALL_O = 0.00001;

    colno[j] = vcol;
    row[j++] = -SMALL_O;

    colno[j] = bcol;
    row[j++] = 1.0;

    add_constraintex(lp, j, row, colno, GE, 0.0);
  }

  for (auto &elem : amounts) {
    std::string const &person = elem.first;
    double amount = elem.second;

    int j = 0;

    if (inflow.find(person) != inflow.end()) {
      for (auto idx : inflow.at(person)) {
        colno[j] = idx;
        row[j++] = 1.0;
      }
    }
    if (outflow.find(person) != outflow.end()) {
      for (auto idx : outflow.at(person)) {
        colno[j] = idx;
        row[j++] = -1.0;
      }
    }
    add_constraintex(lp, j, row, colno, EQ, amount);
  }
  set_add_rowmode(lp, FALSE);

  // Objective functions
  {
    int j = 0;
    for (int idx : binary) {
      colno[j] = idx;
      row[j++] = 1.0;
    }
    set_obj_fnex(lp, j, row, colno);
  }

  set_minim(lp);
  write_LP(lp, stdout);

  set_verbose(lp, IMPORTANT);

  int ret = solve(lp);
  assert(ret == OPTIMAL);

  std::cout << std::endl;
  std::cout << std::endl;
  std::cout << std::endl;

  get_variables(lp, row);
  for (int j : binary) {
    if (row[j - 1] > 0.5) {
      double value = row[j - 2];
      std::pair<std::string, std::string> &payment = payments.at((j - 2) / 2);
      std::cout << payment.first << " to " << payment.second << " " << value
                << std::endl;
    }
  }

  delete_lp(lp);

  return 0;
}
