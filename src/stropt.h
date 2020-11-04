#ifndef STROPT_H
#define STROPT_H

#include <string>
#include <vector>

void    remove_comment(std::string &line);
void    remove_space(std::string &line);                        // remove leading and trailing spaces
void    remove_whitespace(std::string &line);
void    clean_line(std::string &line);

void    str_upper(std::string &line);
void    str_lower(std::string &line);
void    str_change(std::string &line, char pre, char aft);

std::vector<std::string>   str_split(std::string line, std::string flag);



#endif // STROPT_H
