#include <iostream>
#include <map>
#include <vector>
#include <tuple>
#include <algorithm>

/*
 * A helper mapper function for unwrapping an ambiguous base into a vector
 * of its possible values according to IUPAC classification
 * */
std::vector<char> map_bases(char base) {
    std::map<char, std::vector<char>> ambiguous_bases = {
        std::pair<char, std::vector<char>>{'A', std::vector<char>{'A'}},
        std::pair<char, std::vector<char>>{'C', std::vector<char>{'C'}},
        std::pair<char, std::vector<char>>{'G', std::vector<char>{'G'}},
        std::pair<char, std::vector<char>>{'T', std::vector<char>{'T'}},
        std::pair<char, std::vector<char>>{'U', std::vector<char>{'T'}},
        std::pair<char, std::vector<char>>{'N', std::vector<char>{'A', 'C', 'G', 'T'}},
        std::pair<char, std::vector<char>>{'R', std::vector<char>{'A', 'G'}},
        std::pair<char, std::vector<char>>{'Y', std::vector<char>{'C', 'T'}},
        std::pair<char, std::vector<char>>{'S', std::vector<char>{'C', 'G'}},
        std::pair<char, std::vector<char>>{'W', std::vector<char>{'A', 'T'}},
        std::pair<char, std::vector<char>>{'K', std::vector<char>{'G', 'T'}},
        std::pair<char, std::vector<char>>{'M', std::vector<char>{'A', 'C'}},
        std::pair<char, std::vector<char>>{'B', std::vector<char>{'C', 'G', 'T'}},
        std::pair<char, std::vector<char>>{'D', std::vector<char>{'A', 'G', 'T'}},
        std::pair<char, std::vector<char>>{'H', std::vector<char>{'A', 'C', 'T'}},
        std::pair<char, std::vector<char>>{'W', std::vector<char>{'A', 'C', 'G'}},
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
    return (int) intersection.size() > 0 ? match : mismatch;
}

/*
 * Generates a score grid for finding the most optimized resulting sequence
 */
std::vector<std::vector<int>> generate_grid(const std::string &left_sequence,
                                            const std::string &right_sequence,
                                            int match = 1, int mismatch = 0, int gap = -1) {
    std::vector<std::vector<int>> grid;

    // prior assignment of the grid
    std::vector<int> first_row;
    first_row.reserve(left_sequence.length());

    for (int i = 0; i < left_sequence.length(); ++i) {
        first_row.push_back(i * gap);
    }
    grid.push_back(first_row);

    for (int j = 1; j < right_sequence.length(); ++j) {
        std::vector<int> row = {j * gap};
        grid.push_back(row);
    }

    for (int i = 1; i < left_sequence.length(); ++i) {
        for (int j = 1; j < right_sequence.length(); ++j) {
            int cell_match = grid[i - 1][j - 1] + get_score(left_sequence[i], right_sequence[j], match, mismatch);
            int deletion = grid[i - 1][j] + gap;
            int insertion = grid[i][j - 1] + gap;
            grid[i].push_back(std::max({cell_match, deletion, insertion}));
        }
    }

    return grid;
}

std::tuple<std::string, std::string> align(const std::string &seq1,
                                           const std::string &seq2,
                                           int match = 1, int mismatch = 0, int gap = -1) {
    std::string seq_a = "-" + seq1;
    std::string seq_b = "-" + seq2;

    auto grid = generate_grid(seq_a, seq_b, match, mismatch, gap);

    std::string alignment_a, alignment_b;
    int i = (int) seq_a.length() - 1;
    int j = (int) seq_b.length() - 1;

    while (i > 0 && j > 0) {
        int score = grid[i][j];
        int score_diagonal = grid[i - 1][j - 1];
        int score_up = grid[i][j - 1];
        int score_left = grid[i - 1][j];
        if (score == score_diagonal + get_score(seq_a[i], seq_b[j], match, mismatch)) {
            alignment_a += seq_a[i];
            alignment_b += seq_b[j];
            --i;
            --j;
        }
        else if (score == score_left + gap) {
            alignment_a += seq_a[i];
            alignment_b += "-";
            --i;
        }
        else if (score == score_up + gap) {
            alignment_a += "-";
            alignment_b += seq_b[j];
            --j;
        }
    }

    // fill the remain space with gaps
    while (i > 0){
        alignment_a += seq_a[i];
        alignment_b += "-";
        --i;
    }
    while (j > 0){
        alignment_a += "-";
        alignment_b += seq_b[j];
        --j;
    }

    reverse(alignment_a.begin(), alignment_a.end());
    reverse(alignment_b.begin(), alignment_b.end());
    return std::tuple<std::string, std::string>{alignment_a, alignment_b};
}

int main() {
    std::string al1, al2;
    std::tie(al1, al2) = align("GCATGCU", "GATTACA", 1, -1, -1);
    std::cout << al1 << std::endl;
    std::cout << al2 << std::endl;
    return 0;
}
