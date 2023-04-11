#include <iostream>
#include <vector>
#include <iterator>
#include <map>
#include <set>
#include <fstream>
#include <exception>
#include <cassert>

// author @vlburmistrov

namespace IO {
    class Writer {
        int cnt_bits;
        int buffer;
        std::ofstream *outfile;

    public:
        explicit Writer(std::ofstream *file) : cnt_bits(0), buffer(0) {
            outfile = file;
        };
    private:
        void write_in_file(uint8_t x) const {
            outfile->write(reinterpret_cast<const char *>(&x), sizeof x);
        }

    public:
        void clear_buffer() {
            int add = (8 - cnt_bits) % 8;
            if (add == 0) {
                return;
            }
            buffer <<= add;
            write_in_file(buffer);
            cnt_bits = 0;
            buffer = 0;
        }

        void write_bits(std::vector<bool> &num, int n) {
            for (const auto &bit: num) {
                write_bits(bit, 1);
            }
        }

        void write_bits(int num, int n) {
            for (int i = n - 1; i >= 0; i--) {
                int bit = (num & (1 << i)) > 0;
                buffer <<= 1;
                buffer |= bit;
                cnt_bits++;
                if (cnt_bits == 8) {
                    write_in_file(buffer);
                    cnt_bits = 0;
                    buffer = 0;
                }
            }
        }
    };

    /*
    class Reader {
        int ptr;
        int buffer;
        std::ifstream *infile;

    public:
        bool have;

    public:
        explicit Reader(std::ifstream *file) : ptr(7), buffer(0), have(true) {
            infile = file;
        };
    private:
        uint8_t read_byte() {
            char x[1];
            have = (bool) infile->read(x, 1);
            uint8_t num = *reinterpret_cast<uint8_t *>(x);
            return num;
        };
    }
    */
}

class Huffman {
    class Node {
    public:
        unsigned char ch{};
        int cnt{};
        Node *left{};
        Node *right{};

        explicit Node(unsigned char ch, int cnt) : ch(ch), cnt(cnt), left(nullptr), right(nullptr) {}

        Node() = default;

        Node(Node *a, Node *b) {
            left = a;
            right = b;
            cnt = a->cnt + b->cnt;
            ch = std::min(a->ch, b->ch);
        }
    };

    static bool cmp_node(Node *x, Node *y) {
        return x->cnt < y->cnt || (x->cnt == y->cnt && x->ch < y->ch);
    }

    Node *tree{};
    std::map<unsigned char, std::vector<bool>> codes;
    std::vector<unsigned char> data;

public:
    Huffman() = default;

    explicit Huffman(const std::vector<unsigned char> &s) : data(s) {
        build_tree();
        std::vector<bool> vec;
        build_codes(tree, vec);
    }

private:
    void build_codes(Node *root, std::vector<bool> &code) {
        if (root->left == nullptr && root->right == nullptr) {
            codes[root->ch] = (!code.empty() ? code : std::vector<bool>(1));
            return;
        }
        code.push_back(false);
        build_codes(root->left, code);
        code.pop_back();
        code.push_back(true);
        build_codes(root->right, code);
        code.pop_back();
    }

    void build_tree() {
        std::map<unsigned char, int> count;
        for (const unsigned char ch: data) {
            count[ch]++;
        }
        std::multiset<Node *, decltype(&cmp_node)> heap(&cmp_node);
        for (auto &[ch, cnt]: count) {
            heap.insert(new Node(ch, cnt));
        }
        while (heap.size() > 1) {
            auto x = *heap.begin();
            heap.erase(heap.begin());
            auto y = *heap.begin();
            heap.erase(heap.begin());
            heap.insert(new Node(x, y));
        }
        tree = *heap.begin();
    }

    void print_codes() {
        std::cout << codes.size() << std::endl;
        for (auto &[x, y]: codes) {
            std::cout << (int) x << ": ";
            for (bool i: y) {
                std::cout << i;
            }
            std::cout << std::endl;
        }
    }

    void write_tree(Node *root, IO::Writer &w) {
        if (root->left == nullptr && root->right == nullptr) {
            w.write_bits(0, 1);
            w.write_bits(root->ch, 8);
            return;
        }
        w.write_bits(1, 1);
        write_tree(root->left, w);
        write_tree(root->right, w);
    }

public:
    void write_in_file(std::ofstream &out) {
        int cnt_bits = (int) codes.size() * 9 + (int) codes.size() - 1 + 8 * 2 + 3;
        for (const unsigned char ch: data) {
            cnt_bits += (int) codes[ch].size();
        }
        int add = (8 - cnt_bits % 8) % 8;
        cnt_bits += add;
        IO::Writer w(&out);
        w.write_bits('H', 8);
        w.write_bits('A', 8);
        w.write_bits((8 - add) % 8, 3);
        write_tree(tree, w);
        for (const unsigned char ch: data) {
            auto code = codes[ch];
            w.write_bits(code, (int) code.size());
        }
        w.clear_buffer();
    }

    std::vector<unsigned char> decode(std::ifstream &infile) {
        /*
         * тут у меня кончились силы писать нормальный код,
         * поэтому всё будет в одной функции и с моими любимыми лямбдами...
         */
        std::vector<bool> bits;
        bool have;
        auto read_byte = [&have, &infile]() -> uint8_t {
            char x[1];
            have = (bool) infile.read(x, 1);
            uint8_t num = *reinterpret_cast<uint8_t *>(x);
            return num;
        };
        while (true) {
            unsigned char byte = read_byte();
            if (!have) {
                break;
            }
            for (int i = 7; i >= 0; i--) {
                bits.push_back((byte & (1 << i)) > 0);
            }
        }
        int ptr = 0;
        auto bits_to_int = [&ptr, &bits](int len) -> int {
            if (ptr + len - 1 >= (int) bits.size()) {
                throw std::out_of_range("Not enough bits");
            }
            int x = 0;
            for (int it = 0; it < len; it++) {
                x <<= 1;
                x |= bits[ptr++];
            }
            return x;
        };
        unsigned char H = bits_to_int(8);
        unsigned char A = bits_to_int(8);
        assert(H == 72 && A == 65);
        int add = (8 - (bits_to_int(3))) % 8;
        assert((int) bits.size() >= add);
        for (int i = 0; i < add; i++) {
            assert(bits.back() == 0);
            bits.pop_back();
        }
        std::map<std::vector<bool>, unsigned char> char_by_code;
        auto reconstruction_codes = [&bits_to_int, &char_by_code](auto &&reconstruction_codes,
                                                                  std::vector<bool> &code) {
            int bit = bits_to_int(1);
            if (bit == 0) {
                unsigned char ch = bits_to_int(8);
                if (code.empty()) {
                    code = {false};
                }
                char_by_code[code] = ch;
                return;
            }
            code.push_back(false);
            reconstruction_codes(reconstruction_codes, code);
            code.pop_back();
            code.push_back(true);
            reconstruction_codes(reconstruction_codes, code);
            code.pop_back();
        };
        std::vector<bool> code;
        reconstruction_codes(reconstruction_codes, code);
        std::vector<unsigned char> result;
        while (ptr != bits.size()) {
            std::vector<bool> code;
            while (!char_by_code.count(code)) {
                code.push_back(bits_to_int(1));
            }
            result.push_back(char_by_code[code]);
        }
        return result;
    }
};


int main() {
    std::string mode = "encode"; /* "decode" */
    if (mode == "encode") {
        std::ifstream infile("input.dat", std::ios::binary | std::ios::in);
        std::ofstream outfile("output.ha", std::ios::binary | std::ios::out);
        bool have;
        auto read_byte = [&]() -> uint8_t {
            char x[1];
            have = (bool) infile.read(x, 1);
            uint8_t num = *reinterpret_cast<uint8_t *>(x);
            return num;
        };
        std::vector<unsigned char> s;
        while (true) {
            unsigned char ch = read_byte();
            if (!have) {
                break;
            }
            s.push_back(ch);
        }
        Huffman huffman(s);
        huffman.write_in_file(outfile);
    } else {
        std::string input_file_name = "input.ha";
        std::ifstream infile(input_file_name, std::ios::binary | std::ios::in);
        std::ofstream outfile("output.dat", std::ios::binary | std::ios::out);
        Huffman huffman;
        auto data = huffman.decode(infile);
        IO::Writer w(&outfile);
        for (const char ch: data) {
            w.write_bits(ch, 8);
        }
    }
    return 0;
}