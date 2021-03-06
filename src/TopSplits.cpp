#include "TopSplits.h"

void searchGraph(ColoredCDBG<> &graph, unsigned int top_num, string file_name) {

    TopTrie *trie = new TopTrie(top_num);
    ofstream file_out(file_name);
    ostream out(file_out.rdbuf());

    for (auto &unitig : graph) {

        auto num_kmers = unitig.size - Kmer::k + 1;
        auto *uc_kmers = new UnitigColors[num_kmers];

        auto *uc_unitig = unitig.getData()->getUnitigColors(unitig);
        auto unitig_map = UnitigMapBase(0, 1, Kmer::k, true);

        auto it = uc_unitig->begin(unitig);
        auto end = uc_unitig->end();

        for (; it != end; ++it) {
            uc_kmers[it.getKmerPosition()].add(unitig_map, it.getColorID());
        }

        unsigned int len_segment = 1;

        for (unsigned int i = 1; i != num_kmers; ++i, ++len_segment) {

            bool split = (uc_kmers[i].size(unitig_map) != uc_kmers[i - 1].size(unitig_map));

            if (!split) {

                auto curr_it = uc_kmers[i].begin(unitig_map);
                auto curr_end = uc_kmers[i].end();

                auto prev_it = uc_kmers[i - 1].begin(unitig_map);
                auto prev_end = uc_kmers[i - 1].end();

                for (; !split && curr_it != curr_end && prev_it != prev_end; ++curr_it, ++prev_it) {

                    if (curr_it.getColorID() != prev_it.getColorID()) {
                        split = true;
                    }
                }
            }

            if (split) {

                addColors(trie, graph, uc_kmers[i - 1].begin(unitig_map), uc_kmers[i - 1].end(), len_segment);
                len_segment = 0;
            }
        }

        addColors(trie, graph, uc_kmers[num_kmers - 1].begin(unitig_map), uc_kmers[num_kmers - 1].end(), len_segment);
        len_segment = 0;

        delete[] uc_kmers;
    }

    printSplits(trie, out);
}

void addColors(TopTrie *trie, ColoredCDBG<> &graph,
               UnitigColors::const_iterator it, UnitigColors::const_iterator end, unsigned int value) {

    vector<string> graph_colors = graph.getColorNames();
    vector<string> split_colors;

    auto num = distance(it, end);
    auto max = graph_colors.size();

    unsigned int weight = 0;
    unsigned int inverse = 0;

    if (num == 0 || num == max) {
        return;
    } else if (2 * num == max) {
        if (graph_colors.front() == graph.getColorName(it.getColorID())) {
            ++max;
        }
    }

    if (2 * num < max) {
        weight = value;
        for (; it != end; ++it) {
            split_colors.push_back(graph.getColorName(it.getColorID()));
        }
    } else {
        inverse = value;
        string current = graph.getColorName(it.getColorID());
        for (string &color : graph_colors) {
            if (it == end || color != current) {
                split_colors.push_back(color);
            } else {
                current = graph.getColorName((++it).getColorID());
            }
        }
    }

    trie->addNodes(split_colors, weight, inverse);
}

void printSplits(TopTrie *trie, ostream &out) {

    for (auto &split : trie->list) {
        out << split.first;
        for (auto &color : split.second) {
            out << "\t" << color;
        }
        out << endl;
    }
}
