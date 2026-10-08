// Imagine++ project
// Project:  Panorama
// Author:   Pascal Monasse
// Student: TODO

#include <Imagine/Graphics.h>
#include <Imagine/Images.h>
#include <Imagine/LinAlg.h>
#include <vector>
#include <sstream>
using namespace Imagine;
using namespace std;

// Record clicks in two images, until right button click
void getClicks(Window w1, Window w2,
               vector<IntPoint2> &pts1, vector<IntPoint2> &pts2)
{
    // ------------- TODO/A completer ----------
    IntPoint2 p;
    Window win;
    int subWin;
    int pts1_size = 0;
    int pts2_size = 0;

    while (true)
    {
        int button = anyGetMouse(p, win, subWin);

        cout << "Button : " << button << endl;

        if (button == 3)
        {
            if (pts1.size() >= 4 && pts2.size() >= 4)
            {
                cout << "Computing..." << endl
                     << endl;
                break;
            }
            else
            {
                cout << "Please select at least four points for each windows." << endl;
                cout << "Current points window 1: " << pts1.size() << endl;
                cout << "Current points window 2: " << pts2.size() << endl
                     << endl;
            }
        }

        if (win == w1)
        {
            cout << "Window 1" << endl;
            cout << "Coordinate : (" << p.x() << ", " << p.y() << ")" << endl;
            pts1.push_back(p);
            cout << "Number of points selected: " << pts1.size() << endl
                 << endl;
        }
        else if (win == w2)
        {
            cout << "Window 2" << endl;
            cout << "Coordinate : (" << p.x() << ", " << p.y() << ")" << endl;
            pts2.push_back(p);
            cout << "Number of points selected: " << pts2.size() << endl
                 << endl;
        }
    }
}

// Return homography compatible with point matches
Matrix<float> getHomography(const vector<IntPoint2> &pts1,
                            const vector<IntPoint2> &pts2)
{
    size_t n = min(pts1.size(), pts2.size());
    if (n < 4)
    {
        cout << "Not enough correspondences: " << n << endl;
        return Matrix<float>::Identity(3);
    }
    Matrix<double> A(2 * n, 8); // lignes x colonnes
    Vector<double> B(2 * n);
    // ------------- TODO/A completer ----------
    for (size_t i = 0; i < n; i++)
    {
        double x_pts1 = pts1[i].x();
        double y_pts1 = pts1[i].y();

        double x_pts2 = pts2[i].x(); // x'
        double y_pts2 = pts2[i].y(); // y'

        A(2 * i, 0) = x_pts1;
        A(2 * i, 1) = y_pts1;
        A(2 * i, 2) = 1;
        A(2 * i, 3) = 0;
        A(2 * i, 4) = 0;
        A(2 * i, 5) = 0;
        A(2 * i, 6) = -x_pts2 * x_pts1;
        A(2 * i, 7) = -x_pts2 * y_pts1;

        B[2 * i] = x_pts2;

        A(2 * i + 1, 0) = 0;
        A(2 * i + 1, 1) = 0;
        A(2 * i + 1, 2) = 0;
        A(2 * i + 1, 3) = x_pts1;
        A(2 * i + 1, 4) = y_pts1;
        A(2 * i + 1, 5) = 1;
        A(2 * i + 1, 6) = -y_pts2 * x_pts1;
        A(2 * i + 1, 7) = -y_pts2 * y_pts1;

        B[2 * i + 1] = y_pts2;
    }

    B = linSolve(A, B);
    Matrix<float> H(3, 3);
    H(0, 0) = B[0];
    H(0, 1) = B[1];
    H(0, 2) = B[2];
    H(1, 0) = B[3];
    H(1, 1) = B[4];
    H(1, 2) = B[5];
    H(2, 0) = B[6];
    H(2, 1) = B[7];
    H(2, 2) = 1;

    // Sanity check
    for (size_t i = 0; i < n; i++)
    {
        float v1[] = {(float)pts1[i].x(), (float)pts1[i].y(), 1.0f};
        float v2[] = {(float)pts2[i].x(), (float)pts2[i].y(), 1.0f};
        Vector<float> x1(v1, 3);
        Vector<float> x2(v2, 3);
        x1 = H * x1;
        cout << x1[1] * x2[2] - x1[2] * x2[1] << ' '
             << x1[2] * x2[0] - x1[0] * x2[2] << ' '
             << x1[0] * x2[1] - x1[1] * x2[0] << endl;
    }
    return H;
}

// Grow rectangle of corners (x0,y0) and (x1,y1) to include (x,y)
void growTo(float &x0, float &y0, float &x1, float &y1, float x, float y)
{
    if (x < x0)
        x0 = x;
    if (x > x1)
        x1 = x;
    if (y < y0)
        y0 = y;
    if (y > y1)
        y1 = y;
}

// Panorama construction
void panorama(const Image<Color, 2> &I1, const Image<Color, 2> &I2,
              Matrix<float> H)
{
    Vector<float> v(3);
    // on choisit de travailler dans le repère de I2
    float x0 = 0, y0 = 0, x1 = I2.width(), y1 = I2.height();
    for (int i = 0; i < 2; i++)
        for (int j = 0; j < 2; j++)
        {
            v[0] = j * I1.width();
            v[1] = i * I1.height();
            v[2] = 1;
            v = H * v;
            v /= v[2];
            growTo(x0, y0, x1, y1, v[0], v[1]);
        }
    cout << "Rectangle of mosaic in I2 coordinates:" << endl;
    cout << "x0 x1 y0 y1=" << x0 << ' ' << x1 << ' ' << y0 << ' ' << y1 << endl;

    Image<Color> I(int(x1 - x0), int(y1 - y0));
    setActiveWindow(openWindow(I.width(), I.height(), "Panorama"));
    I.fill(WHITE);
    // ------------- TODO/A completer ----------
    Matrix<float> Hinv = inverse(H);

    for (int y = 0; y < I.height(); y++)
    {
        for (int x = 0; x < I.width(); x++)
        {
            float X = x + x0;
            float Y = y + y0;

            // choisir un pixel p2 de l'image I (repère de base: I2)
            Vector<float> p2(3);
            p2[0] = X;
            p2[1] = Y;
            p2[2] = 1;

            /*
                pour obtenir le pixel p1 de I1 depuis le repère de I2
                on applique la formule:
                p1 = H^-1 * p2
                ici on a directement p2, car on utilise P2 comme base de repère
            */
            Vector<float> p1 = Hinv * p2;

            // normaliser
            p1[0] = p1[0] / p1[2];
            p1[1] = p1[1] / p1[2];
            p1[2] = p1[2] / p1[2];

            bool inI1 = p1[0] >= 0 && p1[0] < I1.width() &&
                        p1[1] >= 0 && p1[1] < I1.height();

            bool inI2 = X >= 0 && X < I2.width() &&
                        Y >= 0 && Y < I2.height();

            if (inI1 && inI2)
            {
                // couleur moyenne
                Color c1 = I1(int(p1[0]), int(p1[1]));
                Color c2 = I2(int(X), int(Y));

                Color c(
                    (c1.r() + c2.r()) / 2,
                    (c1.g() + c2.g()) / 2,
                    (c1.b() + c2.b()) / 2);

                I(x, y) = c;
            }
            else if (inI1)
            {
                I(x, y) = I1(int(p1[0]), int(p1[1]));
            }
            else if (inI2)
            {
                I(x, y) = I2(int(X), int(Y));
            }
        }
    }

    display(I, 0, 0);
}

// Main function
int main(int argc, char *argv[])
{
    string s1 = argc > 2 ? argv[1] : srcPath("data/image0006.jpg");
    string s2 = argc > 2 ? argv[2] : srcPath("data/image0007.jpg");

    // Load and display images
    Image<Color> I1, I2;
    if (!load(I1, s1) ||
        !load(I2, s2))
    {
        cerr << "Unable to load the images" << endl;
        return 1;
    }
    Window w1 = openWindow(I1.width(), I1.height(), s1);
    display(I1, 0, 0);
    Window w2 = openWindow(I2.width(), I2.height(), s2);
    setActiveWindow(w2);
    display(I2, 0, 0);

    // Get user's clicks in images
    vector<IntPoint2> pts1, pts2;
    getClicks(w1, w2, pts1, pts2);

    vector<IntPoint2>::const_iterator it;
    cout << "pts1=" << endl;
    for (it = pts1.begin(); it != pts1.end(); it++)
        cout << *it << endl;
    cout << "pts2=" << endl;
    for (it = pts2.begin(); it != pts2.end(); it++)
        cout << *it << endl;

    // Compute homography
    Matrix<float> H = getHomography(pts1, pts2);
    cout << "H=" << H / H(2, 2);

    // Apply homography
    panorama(I1, I2, H);

    endGraphics();
    return 0;
}
