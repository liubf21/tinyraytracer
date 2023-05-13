#pragma once
#include "utils.hpp"

struct PPMnode // photon map node (actually store the intersection point of the ray and the object)
{
    Vec3f pos, col, dir; // position, color, direction
    int index;           // index in the photon map
    float prob, r;       // probability of the photon, radius of the photon
    PPMnode() : index(-1), prob(1){};
    PPMnode(Point3 pos_, Vec3f col_, Vec3f dir_, int index_ = -1, float prob_ = 1, float r_ = 1) : pos(pos_), col(col_), dir(dir_), index(index_), prob(prob_), r(r_){};
};

struct IMGbuf
{
    float n; // number of samples
    Vec3f f; // sum of colors
    IMGbuf() : n(0), f(0, 0, 0) {}
    IMGbuf(float n_, Vec3f f_) : n(n_), f(f_) {}
    void add(Vec3f c, float p = 1.0)
    {
        n += p;
        // std::cout << "n: " << n << std::endl;
        f = f + c * p;
    }
    Vec3f getcol() { return f / n; }
    void clear()
    {
        n = 0;
        f = Vec3f();
    }
    IMGbuf operator+(const IMGbuf &b) const { return IMGbuf(n + b.n, f + b.f); }
    IMGbuf operator*(const float &b) const { return IMGbuf(n * b, f * b); }
    IMGbuf operator/(const float &b) const { return IMGbuf(n / b, f / b); }
    Vec3f get() { return n < epsilon ? f : f / n; }
};

class KDTree
{
public:
    static int D; // dimension used for splitting
    int n, root;
    struct KDTreeNode
    {
        PPMnode pm; // photon map node
        Vec3f m[2]; // min, max
        int s[2];   // left and right son
        KDTreeNode() : pm(), s{-1, -1} {}
        bool operator<(const KDTreeNode &b) const
        {
            return pm.pos[D] < b.pm.pos[D];
        }
    };
    KDTreeNode *tree;
    KDTree() : n(0), root(-1), tree(nullptr) {}
    KDTree(std::vector<PPMnode> &pm) : n(0), root(-1), tree(nullptr) { init(pm); }
    ~KDTree()
    {
        if (tree != nullptr)
            delete[] tree;
    }
    void mt(int f, int x) // modify the min and max of f using x
    {
        // if (f == -1) return;
        tree[f].m[0] = min(tree[f].m[0], tree[x].m[0]);
        tree[f].m[1] = max(tree[f].m[1], tree[x].m[1]);
    }
    int bt(int l, int r, int d)
    {
        D = d;
        int mid = (l + r) >> 1;
        std::nth_element(tree + l, tree + mid, tree + r + 1);
        tree[mid].m[0] = tree[mid].pm.pos - tree[mid].pm.r; // min
        tree[mid].m[1] = tree[mid].pm.pos + tree[mid].pm.r; // max
        if (l < mid)
            tree[mid].s[0] = bt(l, mid - 1, (d + 1) % 3), mt(mid, tree[mid].s[0]); // left
        if (mid < r)
            tree[mid].s[1] = bt(mid + 1, r, (d + 1) % 3), mt(mid, tree[mid].s[1]); // right
        return mid;
    }
    void init(std::vector<PPMnode> &pm) // use PPMnode to initialize the tree
    {
        n = pm.size();
        if (tree != nullptr)
            delete[] tree;
        tree = new KDTreeNode[n + 10]; // make sure there is enough space
        for (int i = 0; i < n; ++i)
            tree[i].pm = pm[i]; // copy
        root = bt(0, n - 1, 0); // build tree
    }
    float getdis2(Vec3f pos, Vec3f m0, Vec3f m1) // get the distance^2 between pos and the box
    {
        return (max(Vec3f(), max(m0 - pos, pos - m1))).norm2();
    }
    void _query(const PPMnode &node, IMGbuf *c, int o) // query whether the node is in the box of o, if so, add to the buffer
    {
        if ((tree[o].pm.pos - node.pos).norm2() <= tree[o].pm.r * tree[o].pm.r) // in the sphere, add to the buffer
        {
            c[tree[o].pm.index].add(wiseProduct(node.col, tree[o].pm.col), node.prob); // calculate the color using the color of the photon
            // std::cout << tree[o].pm.index << "add: " << node.prob <<"->"<< c[tree[o].pm.index].n << std::endl;
        }
        float dis2[2];
        if (tree[o].s[0] != -1)
            dis2[0] = getdis2(node.pos, tree[tree[o].s[0]].m[0], tree[tree[o].s[0]].m[1]);
        else
            dis2[0] = max_float;
        if (tree[o].s[1] != -1)
            dis2[1] = getdis2(node.pos, tree[tree[o].s[1]].m[0], tree[tree[o].s[1]].m[1]);
        else
            dis2[1] = max_float;
        int tmp = dis2[0] >= dis2[1];
        if (dis2[tmp] < epsilon)
            _query(node, c, tree[o].s[tmp]); // if in the box, recursively query
        tmp ^= 1;
        if (dis2[tmp] < epsilon)
            _query(node, c, tree[o].s[tmp]);
    }
    void _modify(int o) // modify the min and max of o, used when the radius of o is changed
    {
        tree[o].m[0] = tree[o].pm.pos - tree[o].pm.r; // min
        tree[o].m[1] = tree[o].pm.pos + tree[o].pm.r; // max
        if (tree[o].s[0] != -1)
            _modify(tree[o].s[0]), mt(o, tree[o].s[0]);
        if (tree[o].s[1] != -1)
            _modify(tree[o].s[1]), mt(o, tree[o].s[1]);
    }
    void query(PPMnode node, IMGbuf *c) // query whether the node is in the box of the root
    {
        _query(node, c, root);
    }
    void modify() // modify the min and max of the root
    {
        _modify(root);
    }
};
