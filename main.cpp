#include <iostream>
#include <map>
#include <vector>
#include <algorithm>

/*
 * A helper mapper function for unwrapping an ambiguous base into a vector
 * of its possible values according to IUPAC classification
 * */
std::vector<char> map_bases(char base) {
    std::map<char, std::vector<char>> ambiguous_bases = {
        std::pair<char, std::vector<char>> {'A', std::vector<char> {'A'}},
        std::pair<char, std::vector<char>> {'C', std::vector<char> {'C'}},
        std::pair<char, std::vector<char>> {'G', std::vector<char> {'G'}},
        std::pair<char, std::vector<char>> {'T', std::vector<char> {'T'}},
        std::pair<char, std::vector<char>> {'N', std::vector<char> {'A', 'C', 'G', 'T'}},
        std::pair<char, std::vector<char>> {'R', std::vector<char> {'A', 'G'}},
        std::pair<char, std::vector<char>> {'Y', std::vector<char> {'C', 'T'}},
        std::pair<char, std::vector<char>> {'S', std::vector<char> {'C', 'G'}},
        std::pair<char, std::vector<char>> {'W', std::vector<char> {'A', 'T'}},
        std::pair<char, std::vector<char>> {'K', std::vector<char> {'G', 'T'}},
        std::pair<char, std::vector<char>> {'M', std::vector<char> {'A', 'C'}},
        std::pair<char, std::vector<char>> {'B', std::vector<char> {'C', 'G', 'T'}},
        std::pair<char, std::vector<char>> {'D', std::vector<char> {'A', 'G', 'T'}},
        std::pair<char, std::vector<char>> {'H', std::vector<char> {'A', 'C', 'T'}},
        std::pair<char, std::vector<char>> {'W', std::vector<char> {'A', 'C', 'G'}},
    };

    return ambiguous_bases.at(base);
}

/*
 * Calculates a matching score for two adjacent sequence elements.
 * If one of the elements is a gap ('.' or '-'), the function will assume a mismatch
 * */
int get_score(char lho, char rho, int match, int mismatch) {
    if (lho == '-' || lho == '.' || rho == '-' || rho == '.') {
        return mismatch;
    }
    auto lho_vector = map_bases(lho);
    auto rho_vector = map_bases(rho);
    std::vector<char> intersection;
    std::set_intersection(lho_vector.begin(), lho_vector.end(),
                          rho_vector.begin(), rho_vector.end(),
                          std::back_inserter(intersection));
    return (int)intersection.size() > 0 ? match : mismatch;
}

int main() {
    std::cout << get_score('A', 'C', 1, 0) << std::endl;
    std::cout << get_score('A', 'N', 1, 0) << std::endl;
    std::cout << get_score('A', 'D', 1, 0) << std::endl;
    std::cout << get_score('A', 'Y', 1, 0) << std::endl;
    return 0;
}
