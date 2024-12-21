// #include <bits/stdc++.h>

#include <cstddef>
#include <cstdint>
#include <iostream>
#include <iomanip>
#include <limits>
#include <array>
#include <fstream>
#include <string>
#include <regex>

using namespace std;

using Size = std::pair<size_t, size_t>;

constexpr size_t N = 36, M = 84;
// constexpr size_t N = 14, M = 5;
constexpr size_t X = 1'000'000;
constexpr std::array<pair<int, int>, 4> deltas{{{-1, 0}, {1, 0}, {0, -1}, {0, 1}}};

// char field[N][M + 1] = {
//     "#####",
//     "#.  #",
//     "#.# #",
//     "#.# #",
//     "#.# #",
//     "#.# #",
//     "#.# #",
//     "#.# #",
//     "#...#",
//     "#####",
//     "#   #",
//     "#   #",
//     "#   #",
//     "#####",
// };

char field[N][M + 1] = {
    "####################################################################################",
    "#                                                                                  #",
    "#                                                                                  #",
    "#                                                                                  #",
    "#                                                                                  #",
    "#                                                                                  #",
    "#                                       .........                                  #",
    "#..............#            #           .........                                  #",
    "#..............#            #           .........                                  #",
    "#..............#            #           .........                                  #",
    "#..............#            #                                                      #",
    "#..............#            #                                                      #",
    "#..............#            #                                                      #",
    "#..............#            #                                                      #",
    "#..............#............#                                                      #",
    "#..............#............#                                                      #",
    "#..............#............#                                                      #",
    "#..............#............#                                                      #",
    "#..............#............#                                                      #",
    "#..............#............#                                                      #",
    "#..............#............#                                                      #",
    "#..............#............#                                                      #",
    "#..............#............################                     #                 #",
    "#...........................#....................................#                 #",
    "#...........................#....................................#                 #",
    "#...........................#....................................#                 #",
    "##################################################################                 #",
    "#                                                                                  #",
    "#                                                                                  #",
    "#                                                                                  #",
    "#                                                                                  #",
    "#                                                                                  #",
    "#                                                                                  #",
    "#                                                                                  #",
    "#                                                                                  #",
    "####################################################################################",
};

template<std::size_t N, std::size_t K, bool F = false>
struct Fixed {
    using T = std::conditional_t<F,
        std::conditional_t<N <= 8, std::int_fast8_t,
        std::conditional_t<N <= 16, std::int_fast16_t,
        std::conditional_t<N <= 32, std::int_fast32_t,
        std::conditional_t<N <= 64, std::int_fast64_t, void>>>>,
        std::conditional_t<N == 8, std::int8_t,
        std::conditional_t<N == 16, std::int16_t, 
        std::conditional_t<N == 32, std::int32_t,
        std::conditional_t<N == 64, std::int64_t, void>>>>>;
    
    static constexpr std::size_t n = N;
    static constexpr std::size_t k = K;
    static constexpr bool f = F;

    constexpr Fixed(): v(0) {}
    constexpr Fixed(float o) {
        v = static_cast<T>(o * (1LL << K));
    }
    constexpr Fixed(double o) {
        v = static_cast<T>(o * (1LL << K));
    }
    template<std::size_t N_, std::size_t K_, bool F_>
    constexpr Fixed(Fixed<N_, K_, F_> o) {
        if (K >= K_) v = static_cast<T>(o.v) * (static_cast<T>(1) << (K - K_));
        else v = static_cast<T>(o.v) / (static_cast<T>(1) << (K_ - K));
    }

    static constexpr Fixed from_raw(T x) {
        Fixed ret;
        ret.v = x;
        return ret;
    }

    operator float() const {
        return static_cast<float>(v) / (1LL << K);
    }

    operator double() const {
        return static_cast<double>(v) / (1LL << K);
    }

    T v;

    auto operator<=>(const Fixed&) const = default;
    bool operator==(const Fixed&) const = default;
};

template<std::size_t N, std::size_t K>
using FastFixed = Fixed<N, K, true>;

template<std::size_t N, std::size_t K, bool F>
static constexpr Fixed<N, K, F> inf = Fixed<N, K, F>::from_raw(std::numeric_limits<typename Fixed<N, K, F>::T>::max());
// template<std::size_t N, std::size_t K, bool F>
// static constexpr Fixed<N, K, F> eps = Fixed<N, K, F>::from_raw(deltas.size());

template<std::size_t N, std::size_t K, bool F>
Fixed<N, K, F> operator+(Fixed<N, K, F> a, Fixed<N, K, F> b) {
    return Fixed<N, K, F>::from_raw(a.v + b.v);
}

template<std::size_t N, std::size_t K, bool F>
Fixed<N, K, F> operator-(Fixed<N, K, F> a, Fixed<N, K, F> b) {
    return Fixed<N, K, F>::from_raw(a.v - b.v);
}

template<std::size_t N, std::size_t K, bool F>
Fixed<N, K, F> operator*(Fixed<N, K, F> a, Fixed<N, K, F> b) {
    return Fixed<N, K, F>::from_raw(((int64_t) a.v * b.v) >> K);
}

template<std::size_t N, std::size_t K, bool F>
Fixed<N, K, F> operator/(Fixed<N, K, F> a, Fixed<N, K, F> b) {
    return Fixed<N, K, F>::from_raw(((int64_t) a.v << K) / b.v);
}

template<std::size_t N, std::size_t K, bool F>
Fixed<N, K, F> &operator+=(Fixed<N, K, F> &a, Fixed<N, K, F> b) {
    return a = a + b;
}

template<std::size_t N, std::size_t K, bool F>
Fixed<N, K, F> &operator-=(Fixed<N, K, F> &a, Fixed<N, K, F> b) {
    return a = a - b;
}

template<std::size_t N, std::size_t K, bool F>
Fixed<N, K, F> &operator*=(Fixed<N, K, F> &a, Fixed<N, K, F> b) {
    return a = a * b;
}

template<std::size_t N, std::size_t K, bool F>
Fixed<N, K, F> &operator/=(Fixed<N, K, F> &a, Fixed<N, K, F> b) {
    return a = a / b;
}

template<std::size_t N, std::size_t K, bool F>
Fixed<N, K, F> operator-(Fixed<N, K, F> x) {
    return Fixed<N, K, F>::from_raw(-x.v);
}

template<std::size_t N, std::size_t K, bool F>
Fixed<N, K, F> abs(Fixed<N, K, F> x) {
    if (x.v < 0) -x;
    return x;
}

template<std::size_t N, std::size_t K, bool F>
std::ostream &operator<<(std::ostream &out, Fixed<N, K, F> x) {
    return out << (double)x;
}

template<std::size_t N, std::size_t M, typename T>
struct VectorField {
    std::array<T, 4> v[N][M];
    T &add(int x, int y, int dx, int dy, T dv) {
        return get(x, y, dx, dy) += dv;
    }
    T &get(int x, int y, int dx, int dy) {
        std::size_t i = std::ranges::find(deltas, std::pair(dx, dy)) - deltas.begin();
        assert(i < deltas.size());
        return v[x][y][i];
    }
};

template<typename T>
struct ParticleParams {
    char type;
    T cur_p;
    array<T, deltas.size()> v;

    void swap_with(int x, int y) {
        swap(field[x][y], type);
        swap(p[x][y], cur_p);
        swap(velocity.v[x][y], v);
    }
};

// template<std::size_t N, std::size_t M, typename T1, typename T2, typename T3>
// struct State {
    
// };

template<std::size_t N, std::size_t M, typename T1, typename T2, typename T3>
class Simulator {
private:
    T1 rho[256];
    T1 p[N][M]{};
    T1 old_p[N][M]{};

    VectorField<N, M, T2> veclocity{}, velocity_flow{};

    int last_use[N][M]{};
    int UT = 0;

    std::mt19937 rnd(1337);

    int dirs[N][M]{};

    tuple<T2, bool, pair<int, int>> propagate_flow(int x, int y, T2 lim) {
        last_use[x][y] = UT - 1;
        T2 ret = 0;
        for (auto [dx, dy] : deltas) {
            int nx = x + dx, ny = y + dy;
            if (field[nx][ny] != '#' && last_use[nx][ny] < UT) {
                auto cap = velocity.get(x, y, dx, dy);
                auto flow = velocity_flow.get(x, y, dx, dy);
                if (flow == cap) {
                    continue;
                }
                
                auto vp = min(lim, cap - flow);
                if (last_use[nx][ny] == UT - 1) {
                    velocity_flow.add(x, y, dx, dy, vp);
                    last_use[x][y] = UT;
                    return {vp, 1, {nx, ny}};
                }
                auto [t, prop, end] = propagate_flow(nx, ny, vp);
                ret += t;
                if (prop) {
                    velocity_flow.add(x, y, dx, dy, t);
                    last_use[x][y] = UT;
                    return {t, prop && end != std::pair(x, y), end};
                }
            }
        }
        last_use[x][y] = UT;
        return {ret, 0, {0, 0}};
    }

    T2 random01() {
        return T2::from_raw((rnd() & ((1 << (T2::n / 2)) - 1)));
    }

    void propagate_stop(int x, int y, bool force = false) {
        if (!force) {
            bool stop = true;
            for (auto [dx, dy] : deltas) {
                int nx = x + dx, ny = y + dy;
                if (field[nx][ny] != '#' && last_use[nx][ny] < UT - 1 && velocity.get(x, y, dx, dy) > 0) {
                    stop = false;
                    break;
                }
            }
            if (!stop) {
                return;
            }
        }
        last_use[x][y] = UT;
        for (auto [dx, dy] : deltas) {
            int nx = x + dx, ny = y + dy;
            if (field[nx][ny] == '#' || last_use[nx][ny] == UT || velocity.get(x, y, dx, dy) > 0) {
                continue;
            }
            propagate_stop(nx, ny);
        }
    }

    T2 move_prob(int x, int y) {
        T2 sum = 0;
        for (size_t i = 0; i < deltas.size(); ++i) {
            auto [dx, dy] = deltas[i];
            int nx = x + dx, ny = y + dy;
            if (field[nx][ny] == '#' || last_use[nx][ny] == UT) {
                continue;
            }
            auto v = velocity.get(x, y, dx, dy);
            if (v < 0) {
                continue;
            }
            sum += v;
        }
        return sum;
    }

    bool propagate_move(int x, int y, bool is_first) {
        last_use[x][y] = UT - is_first;
        bool ret = false;
        int nx = -1, ny = -1;
        do {
            std::array<Fixed, deltas.size()> tres;
            T2 sum = 0;
            for (size_t i = 0; i < deltas.size(); ++i) {
                auto [dx, dy] = deltas[i];
                int nx = x + dx, ny = y + dy;
                if (field[nx][ny] == '#' || last_use[nx][ny] == UT) {
                    tres[i] = sum;
                    continue;
                }
                auto v = velocity.get(x, y, dx, dy);
                if (v < 0) {
                    tres[i] = sum;
                    continue;
                }
                sum += v;
                tres[i] = sum;
            }

            if (sum == (T2)0.) {
                break;
            }

            T2 p = random01() * sum;
            size_t d = std::ranges::upper_bound(tres, p) - tres.begin();

            auto [dx, dy] = deltas[d];
            nx = x + dx;
            ny = y + dy;
            assert(velocity.get(x, y, dx, dy) > (T2)0. && field[nx][ny] != '#' && last_use[nx][ny] < UT);

            ret = (last_use[nx][ny] == UT - 1 || propagate_move(nx, ny, false));
        } while (!ret);
        last_use[x][y] = UT;
        for (size_t i = 0; i < deltas.size(); ++i) {
            auto [dx, dy] = deltas[i];
            int nx = x + dx, ny = y + dy;
            if (field[nx][ny] != '#' && last_use[nx][ny] < UT - 1 && velocity.get(x, y, dx, dy) < (T2)0.) {
                propagate_stop(nx, ny);
            }
        }
        if (ret) {
            if (!is_first) {
                ParticleParams<T1> pp{};
                pp.swap_with(x, y);
                pp.swap_with(nx, ny);
                pp.swap_with(x, y);
            }
        }
        return ret;
    }
public:
    Simulator() {}

    void simulate() {
        rho[' '] = (T1)0.01;
        rho['.'] = (T1)1000.;

        T1 g = (T1)0.1;

        for (size_t x = 0; x < N; ++x) {
            for (size_t y = 0; y < M; ++y) {
                if (field[x][y] == '#')
                    continue;
                for (auto [dx, dy] : deltas) {
                    dirs[x][y] += (field[x + dx][y + dy] != '#');
                }
            }
        }

        for (size_t i = 0; i < X; ++i) {
            T1 total_delta_p = 0.;

            for (size_t x = 0; x < N; ++x) {
                for (size_t y = 0; y < M; ++y) {
                    if (field[x][y] == '#')
                        continue;
                    if (field[x + 1][y] != '#')
                        velocity.add(x, y, 1, 0, g);
                }
            }

            memcpy(old_p, p, sizeof(p));
            for (size_t x = 0; x < N; ++x) {
                for (size_t y = 0; y < M; ++y) {
                    if (field[x][y] == '#')
                        continue;
                    for (auto [dx, dy] : deltas) {
                        int nx = x + dx, ny = y + dy;
                        if (field[nx][ny] != '#' && old_p[nx][ny] < old_p[x][y]) {
                            auto delta_p = old_p[x][y] - old_p[nx][ny];
                            auto force = (T1)delta_p;
                            auto &contr = velocity.get(nx, ny, -dx, -dy);
                            if ((T1)contr * rho[(int) field[nx][ny]] >= force) {
                                contr -= (T2)(force / rho[(int) field[nx][ny]]);
                                continue;
                            }
                            force -= (T1)contr * rho[(int) field[nx][ny]];
                            contr = 0;
                            velocity.add(x, y, dx, dy, (T2)(force / rho[(int) field[x][y]]));
                            p[x][y] -= force / dirs[x][y];
                            total_delta_p -= force / dirs[x][y];
                        }
                    }
                }
            }

            velocity_flow = {};
            bool prop = false;
            do {
                UT += 2;
                prop = 0;
                for (size_t x = 0; x < N; ++x) {
                    for (size_t y = 0; y < M; ++y) {
                        if (field[x][y] != '#' && last_use[x][y] != UT) {
                            auto [t, local_prop, _] = propagate_flow(x, y, 1);
                            if (t > 0) {
                                prop = 1;
                            }
                        }
                    }
                }
            } while (prop);

            for (size_t x = 0; x < N; ++x) {
                for (size_t y = 0; y < M; ++y) {
                    if (field[x][y] == '#')
                        continue;
                    for (auto [dx, dy] : deltas) {
                        auto old_v = velocity.get(x, y, dx, dy);
                        auto new_v = velocity_flow.get(x, y, dx, dy);
                        if (old_v > 0) {
                            assert(new_v <= old_v);
                            velocity.get(x, y, dx, dy) = new_v;
                            auto force = (T1)(old_v - new_v) * rho[(int) field[x][y]];
                            if (field[x][y] == '.')
                                force *= 0.8;
                            if (field[x + dx][y + dy] == '#') {
                                p[x][y] += force / dirs[x][y];
                                total_delta_p += force / dirs[x][y];
                            } else {
                                p[x + dx][y + dy] += force / dirs[x + dx][y + dy];
                                total_delta_p += force / dirs[x + dx][y + dy];
                            }
                        }
                    }
                }
            }

            UT += 2;
            prop = false;
            for (size_t x = 0; x < N; ++x) {
                for (size_t y = 0; y < M; ++y) {
                    if (field[x][y] != '#' && last_use[x][y] != UT) {
                        if (random01() < move_prob(x, y)) {
                            prop = true;
                            propagate_move(x, y, true);
                        } else {
                            propagate_stop(x, y, true);
                        }
                    }
                }
            }

            if (prop) {
                cout << "Tick " << i << ":\n";
                for (size_t x = 0; x < N; ++x) {
                    cout << field[x] << "\n";
                }
            }
        }
    }
};


std::string get_arg(std::string_view arg_name, int argc, char** argv) {
    for (int i = 1; i < argc - 1; ++i) {
        if (argv[i] == arg_name) {
            return argv[i + 1];
        }
    }

    std::cerr << "Error: " << arg_name << " not found." << std::endl;
    exit(1);
}


std::string read_file(const std::string& filename, size_t &n, size_t &m) {
    std::ifstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Error: can't open " << filename << std::endl;
        exit(1);
    }
    
    n = 0; m = 0;

    std::string result, line;
    while (std::getline(file, line)) {
        result += line + '\n';
        n++;
        if (m == 0) m = line.length();
        else if (m == line.length()) continue;

        std::cerr << "Error: lines in file have multiple length" << std::endl;
        exit(1);
    }
    return result;
}


bool check_type(const std::string& type) {
    const std::regex float_regex(R"(^FLOAT$)");
    const std::regex double_regex(R"(^DOUBLE$)");
    const std::regex fixed_regex(R"(^FIXED\((\d+),\s*(\d+)\)$)");
    const std::regex fast_fixed_regex(R"(^FAST_FIXED\((\d+),\s*(\d+)\)$)");

    if (std::regex_match(type, float_regex)) return true;
    if (std::regex_match(type, double_regex)) return true;
    if (std::regex_match(type, fixed_regex)) {
        std::smatch match;
        if (std::regex_match(type, match, fixed_regex)) {
            try {
                int N = std::stoi(match[1].str());
                int K = std::stoi(match[2].str());
                if (N > 0 && K >= 0 && K <= N) return true;
                else return false;
            } catch (const std::invalid_argument&) {
                return false;
            }
        }
    }
    if (std::regex_match(type, fast_fixed_regex)) {
        std::smatch match;
        if (std::regex_match(type, match, fast_fixed_regex)) {
            try {
                int N = std::stoi(match[1].str());
                int K = std::stoi(match[2].str());
                
                if (N > 0 && K >= 0) return true;
                else return false;
            } catch (const std::invalid_argument&) {
                return false;
            }
        }
    }
    
    return false;
}


int main(int argc, char** argv) {
    std::string t1 = get_arg("--p-type", argc, argv);
    std::string t2 = get_arg("--v-type", argc, argv);
    std::string t3 = get_arg("--v-flow-type", argc, argv);
    std::string filename = get_arg("--field", argc, argv);

    if (!check_type(t1)) {
        cerr << "Error: invalid p_type" << endl;
        return 1;
    }
    if (!check_type(t2)) {
        cerr << "Error: invalid v_type" << endl;
        return 1;
    }
    if (!check_type(t3)) {
        cerr << "Error: invalid v-flow_type" << endl;
        return 1;
    }

    size_t n, m;
    std::string field = read_file(filename, n, m);

    
}
