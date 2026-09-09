#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <algorithm>
#include <iomanip>
#include <cctype>

struct ScoringScheme {
    int match = 1;
    int mismatch = -1;
    int gap = -2;
};

struct AlignmentResult {
    int score;
    std::string aligned1;
    std::string aligned2;
};

struct SequenceEntry {
    std::string name;
    std::string seq;
};

AlignmentResult needleman_wunsch(const std::string& seq1, 
                                 const std::string& seq2, 
                                 const ScoringScheme& sc) {
    const size_t n = seq1.size();
    const size_t m = seq2.size();

    std::vector<std::vector<int>> dp(n + 1, std::vector<int>(m + 1, 0));

    for (size_t i = 0; i <= n; ++i) dp[i][0] = static_cast<int>(i) * sc.gap;
    for (size_t j = 0; j <= m; ++j) dp[0][j] = static_cast<int>(j) * sc.gap;

    for (size_t i = 1; i <= n; ++i) {
        for (size_t j = 1; j <= m; ++j) {
            int diag = dp[i - 1][j - 1] + (seq1[i - 1] == seq2[j - 1] ? sc.match : sc.mismatch);
            int up   = dp[i - 1][j] + sc.gap;
            int left = dp[i][j - 1] + sc.gap;
            dp[i][j] = std::max({diag, up, left});
        }
    }

    std::string al1 = "";
    std::string al2 = "";
    size_t i = n;
    size_t j = m;

    while (i > 0 || j > 0) {
        if (i > 0 && dp[i][j] == dp[i - 1][j] + sc.gap) {
            al1 += seq1[--i];
            al2 += '-';
        } else if (j > 0 && dp[i][j] == dp[i][j - 1] + sc.gap) {
            al1 += '-';
            al2 += seq2[--j];
        } else if (i > 0 && j > 0) {
            int match_mismatch = (seq1[i - 1] == seq2[j - 1]) ? sc.match : sc.mismatch;
            if (dp[i][j] == dp[i - 1][j - 1] + match_mismatch) {
                al1 += seq1[--i];
                al2 += seq2[--j];
            }
        }
    }

    std::reverse(al1.begin(), al1.end());
    std::reverse(al2.begin(), al2.end());
    return {dp[n][m], al1, al2};
}

int compute_sp_score(const std::vector<std::string>& msa, const ScoringScheme& sc) {
    int total_score = 0;
    const size_t k = msa.size();
    const size_t len = msa[0].size();

    for (size_t col = 0; col < len; ++col) {
        for (size_t i = 0; i < k; ++i) {
            for (size_t j = i + 1; j < k; ++j) {
                char a = msa[i][col];
                char b = msa[j][col];
                if (a == '-' && b == '-') {
                    continue;
                } else if (a == '-' || b == '-') {
                    total_score += sc.gap;
                } else if (a == b) {
                    total_score += sc.match;
                } else {
                    total_score += sc.mismatch;
                }
            }
        }
    }
    return total_score;
}

struct StarMSAResult {
    size_t center_idx;
    std::vector<std::string> msa;
    int sp_score;
    std::vector<std::vector<int>> pair_scores;
    std::vector<int> row_sums;
};

StarMSAResult star_msa(const std::vector<std::string>& sequences, const ScoringScheme& sc) {
    const size_t k = sequences.size();
    std::vector<std::vector<int>> pair_scores(k, std::vector<int>(k, 0));
    std::vector<int> row_sums(k, 0);

    for (size_t i = 0; i < k; ++i) {
        for (size_t j = i + 1; j < k; ++j) {
            AlignmentResult res = needleman_wunsch(sequences[i], sequences[j], sc);
            pair_scores[i][j] = res.score;
            pair_scores[j][i] = res.score;
        }
    }

    for (size_t i = 0; i < k; ++i) {
        for (size_t j = 0; j < k; ++j) {
            if (i != j) row_sums[i] += pair_scores[i][j];
        }
    }

    size_t center_idx = 0;
    int max_sum = row_sums[0];
    for (size_t i = 1; i < k; ++i) {
        if (row_sums[i] > max_sum) {
            max_sum = row_sums[i];
            center_idx = i;
        }
    }

    std::vector<std::string> current_msa(k, "");
    current_msa[center_idx] = sequences[center_idx];
    bool first_alignment = true;

    for (size_t j = 0; j < k; ++j) {
        if (j == center_idx) continue;

        AlignmentResult pw = needleman_wunsch(sequences[center_idx], sequences[j], sc);
        const std::string& cand_c = pw.aligned1;
        const std::string& cand_s = pw.aligned2;

        if (first_alignment) {
            current_msa[center_idx] = cand_c;
            current_msa[j] = cand_s;
            first_alignment = false;
        } else {
            std::vector<std::string> updated_msa(k, "");
            std::string updated_new_s = "";

            size_t pA = 0;
            size_t pB = 0;

            while (pA < current_msa[center_idx].size() || pB < cand_c.size()) {
                if (pA < current_msa[center_idx].size() && pB < cand_c.size()) {
                    char cA = current_msa[center_idx][pA];
                    char cB = cand_c[pB];

                    if (cA == '-') {
                        for (size_t r = 0; r < k; ++r) {
                            if (!current_msa[r].empty()) updated_msa[r] += current_msa[r][pA];
                        }
                        updated_new_s += '-';
                        pA++;
                    } else if (cB == '-') {
                        for (size_t r = 0; r < k; ++r) {
                            if (!current_msa[r].empty()) updated_msa[r] += '-';
                        }
                        updated_new_s += cand_s[pB];
                        pB++;
                    } else {
                        for (size_t r = 0; r < k; ++r) {
                            if (!current_msa[r].empty()) updated_msa[r] += current_msa[r][pA];
                        }
                        updated_new_s += cand_s[pB];
                        pA++;
                        pB++;
                    }
                } else if (pA < current_msa[center_idx].size()) {
                    for (size_t r = 0; r < k; ++r) {
                        if (!current_msa[r].empty()) updated_msa[r] += current_msa[r][pA];
                    }
                    updated_new_s += '-';
                    pA++;
                } else {
                    for (size_t r = 0; r < k; ++r) {
                        if (!current_msa[r].empty()) updated_msa[r] += '-';
                    }
                    updated_new_s += cand_s[pB];
                    pB++;
                }
            }
            updated_msa[j] = updated_new_s;
            current_msa = updated_msa;
        }
    }

    int sp = compute_sp_score(current_msa, sc);
    return {center_idx, current_msa, sp, pair_scores, row_sums};
}

std::vector<SequenceEntry> load_sequences_from_file(const std::string& filepath) {
    std::ifstream file(filepath);
    if (!file.is_open()) {
        std::cerr << "Error: Cannot open file '" << filepath << "'\n";
        exit(1);
    }

    std::vector<SequenceEntry> entries;
    std::string line;

    while (std::getline(file, line)) {
        if (line.empty()) continue;
        std::istringstream iss(line);
        std::string name, raw_seq;
        if (iss >> name >> raw_seq) {
            std::string clean_seq = "";
            for (char ch : raw_seq) {
                if (std::isalpha(static_cast<unsigned char>(ch))) {
                    clean_seq += std::toupper(static_cast<unsigned char>(ch));
                }
            }
            if (!clean_seq.empty()) {
                entries.push_back({name, clean_seq});
            }
        }
    }

    if (entries.empty()) {
        std::cerr << "Error: No valid sequences found in '" << filepath << "'\n";
        exit(1);
    }

    return entries;
}

int main(int argc, char* argv[]) {
    if (argc < 2) {
        std::cout << "Usage: " << argv[0] << " <input_file.txt>\n";
        return 1;
    }

    std::string input_file = argv[1];
    std::vector<SequenceEntry> entries = load_sequences_from_file(input_file);

    std::vector<std::string> sequences;
    std::vector<std::string> names;
    for (const auto& entry : entries) {
        names.push_back(entry.name);
        sequences.push_back(entry.seq);
    }

    ScoringScheme sc = {1, -1, -2};
    const size_t k = sequences.size();

    StarMSAResult result = star_msa(sequences, sc);

    std::cout << "===============================================================\n";
    std::cout << " 1. PAIRWISE SIMILARITY MATRIX (Needleman-Wunsch)\n";
    std::cout << "===============================================================\n\n";

    std::cout << std::setw(8) << " ";
    for (size_t i = 0; i < k; ++i) std::cout << std::setw(7) << names[i];
    std::cout << std::setw(10) << "Sum\n";

    for (size_t i = 0; i < k; ++i) {
        std::cout << std::setw(8) << names[i];
        for (size_t j = 0; j < k; ++j) {
            std::cout << std::setw(7) << result.pair_scores[i][j];
        }
        std::cout << std::setw(10) << result.row_sums[i] << "\n";
    }

    std::cout << "\nCenter Sequence Selected: " << names[result.center_idx]
              << " with cumulative score = " << result.row_sums[result.center_idx] << "\n\n";

    std::cout << "===============================================================\n";
    std::cout << " 2. RESULTING MULTIPLE SEQUENCE ALIGNMENT (MSA)\n";
    std::cout << "===============================================================\n\n";

    for (size_t i = 0; i < k; ++i) {
        std::cout << std::setw(6) << names[i] << ": ";
        for (char c : result.msa[i]) {
            std::cout << c << "  ";
        }
        std::cout << "\n";
    }

    std::cout << "\n---------------------------------------------------------------\n";
    std::cout << "Final Sum-of-Pairs Score (SP-Score): " << result.sp_score << "\n";
    std::cout << "Alignment Length: " << result.msa[0].size() << " columns\n";
    std::cout << "---------------------------------------------------------------\n";

    return 0;
}