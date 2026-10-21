// Imagine++ project
// Project:  Fundamental
// Author:   Pascal Monasse

#include "./Imagine/Features.h"
#include <Imagine/Graphics.h>
#include <Imagine/LinAlg.h>
#include <vector>
#include <cstdlib>
#include <ctime>
#include <algorithm>
#include <iostream>
#include <random>
using namespace Imagine;
using namespace std;

static const float BETA = 0.01f; // Probability of failure

struct Match
{
    float x1, y1, x2, y2;
};

// Display SIFT points and fill vector of point correspondences
void algoSIFT(Image<Color, 2> I1, Image<Color, 2> I2,
              vector<Match> &matches)
{
    // Find interest points
    SIFTDetector D;
    D.setFirstOctave(-1);
    Array<SIFTDetector::Feature> feats1 = D.run(I1);
    drawFeatures(feats1, Coords<2>(0, 0));
    cout << "Im1: " << feats1.size() << flush;
    Array<SIFTDetector::Feature> feats2 = D.run(I2);
    drawFeatures(feats2, Coords<2>(I1.width(), 0));
    cout << " Im2: " << feats2.size() << flush;

    const double MAX_DISTANCE = 100.0 * 100.0;
    for (size_t i = 0; i < feats1.size(); i++)
    {
        SIFTDetector::Feature f1 = feats1[i];
        for (size_t j = 0; j < feats2.size(); j++)
        {
            double d = squaredDist(f1.desc, feats2[j].desc);
            if (d < MAX_DISTANCE)
            {
                Match m;
                m.x1 = f1.pos.x();
                m.y1 = f1.pos.y();
                m.x2 = feats2[j].pos.x();
                m.y2 = feats2[j].pos.y();
                matches.push_back(m);
            }
        }
    }
}

// RANSAC algorithm to compute F from point matches (8-point algorithm)
// Parameter matches is filtered to keep only inliers as output.
FMatrix<float, 3, 3> computeF(vector<Match> &matches)
{
    const float distMax = 1.5f; // Pixel error for inlier/outlier discrimination
    int Niter = 100000;         // Adjusted dynamically
    FMatrix<float, 3, 3> bestF;
    vector<int> bestInliers;
    // --------------- TODO ------------
    // DO NOT FORGET NORMALIZATION OF POINTS
    float NORM_util[3][3] =
        {
            {0.001, 0, 0},
            {0, 0.001, 0},
            {0, 0, 1}};

    FMatrix<float, 3, 3> NORM(NORM_util);

    // 1. normaliser les coordonnées
    vector<Match> normalized_matches;
    for (int i = 0; i < matches.size(); i++)
    {
        Match m = matches[i];
        m.x1 = m.x1 * 0.001;
        m.y1 = m.y1 * 0.001;
        m.x2 = m.x2 * 0.001;
        m.y2 = m.y2 * 0.001;

        normalized_matches.push_back(m);
    }

    mt19937 gen(random_device{}());

    // RANSAC
    for (int iter = 0; iter < Niter; iter++)
    {
        // 2. sélectionner 8 matches différents
        vector<Match> out;
        sample(
            normalized_matches.begin(),
            normalized_matches.end(),
            back_inserter(out),
            8,
            gen);

        // 3. construire A
        FMatrix<float, 9, 9> A;

        for (int i = 0; i < 8; i++)
        {
            Match m = out[i];

            A(i, 0) = m.x1 * m.x2;
            A(i, 1) = m.x1 * m.y2;
            A(i, 2) = m.x1;
            A(i, 3) = m.y1 * m.x2;
            A(i, 4) = m.y1 * m.y2;
            A(i, 5) = m.y1;
            A(i, 6) = m.x2;
            A(i, 7) = m.y2;
            A(i, 8) = 1;
        }

        // remplir la dernière ligne avec des zéros pour une matrice de 9x9
        for (int j = 0; j < 9; j++)
        {
            A(8, j) = 0;
        }

        // 4. SVD(A) -> F
        FMatrix<float, 9, 9> U;
        FVector<float, 9> S;
        FMatrix<float, 9, 9> Vt;
        svd(A, U, S, Vt);

        FMatrix<float, 3, 3> F;

        for (int i = 0; i < 3; i++)
        {
            for (int j = 0; j < 3; j++)
            {
                F(i, j) = Vt(8, 3 * i + j);
            }
        }

        // 5. imposer rank(F)=2
        FMatrix<float, 3, 3> Uf;
        FVector<float, 3> Sf;
        FMatrix<float, 3, 3> Vtf;

        svd(F, Uf, Sf, Vtf);

        Sf[2] = 0;

        FMatrix<float, 3, 3> D;
        // ligne 1
        D(0, 0) = Sf[0];
        D(0, 1) = 0;
        D(0, 2) = 0;
        // ligne 2
        D(1, 0) = 0;
        D(1, 1) = Sf[1];
        D(1, 2) = 0;
        // ligne 2
        D(2, 0) = 0;
        D(2, 1) = 0;
        D(2, 2) = 0;

        F = Uf * D * Vtf;

        // 6. dénormaliser F
        F = transpose(NORM) * F * NORM;

        // 7. tester tous les matches
        vector<int> inliers;

        for (int i = 0; i < matches.size(); i++)
        {
            Match m = matches[i];

            float x = m.x1;
            float y = m.y1;
            float xp = m.x2;
            float yp = m.y2;

            // droite épipolaire dans l'image 2 : l = F^T * x
            float a = F(0, 0) * x + F(1, 0) * y + F(2, 0);
            float b = F(0, 1) * x + F(1, 1) * y + F(2, 1);
            float c = F(0, 2) * x + F(1, 2) * y + F(2, 2);

            float denom = sqrt(a * a + b * b);

            if (denom > 0)
            {
                float d = abs(a * xp + b * yp + c) / denom;

                if (d < distMax)
                {
                    inliers.push_back(i);
                }
            }
        }

        // 8. conserver les inliers si meilleur modèle
        if (inliers.size() > bestInliers.size())
        {
            bestInliers = inliers;
            bestF = F;

            // 9. mettre Niter à jour
            float ratio = float(bestInliers.size()) / float(matches.size());

            float pBadSample = 1.0f - pow(ratio, 8);

            if (pBadSample > 0.0f && pBadSample < 1.0f)
            {
                int newNiter =
                    int(ceil(log(BETA) / log(pBadSample)));

                if (newNiter < Niter)
                {
                    Niter = newNiter;
                }
            }
        }
    }

    // Updating matches with inliers only
    vector<Match> all = matches;
    matches.clear();
    for (size_t i = 0; i < bestInliers.size(); i++)
        matches.push_back(all[bestInliers[i]]);
    return bestF;
}

// Expects clicks in one image and show corresponding line in other image.
// Stop at right-click.
void displayEpipolar(Image<Color> I1, Image<Color> I2,
                     const FMatrix<float, 3, 3> &F)
{
    int w = I1.width();
    int h = I1.height();

    while (true)
    {
        int x, y;

        if (getMouse(x, y) == 3)
            break;

        Color c(rand() % 256, rand() % 256, rand() % 256);
        fillCircle(x, y, 2, c);

        float p_util[3][1] = {{(float)x}, {(float)y}, {1}};
        FMatrix<float, 3, 1> p(p_util);

        FMatrix<float, 3, 1> l;
        int decalage = 0;

        if (x < w)
        {
            l = transpose(F) * p;
            decalage = w;
        }
        else
        {
            p(0, 0) -= w;
            l = F * p;
        }

        float a = l(0, 0);
        float b = l(1, 0);
        float cc = l(2, 0);

        vector<float> pts;

        if (b != 0)
        {
            float y1 = -cc / b;
            if (y1 >= 0 && y1 <= h)
            {
                pts.push_back(decalage);
                pts.push_back(y1);
            }

            float y2 = (-cc - a * w) / b;
            if (y2 >= 0 && y2 <= h)
            {
                pts.push_back(w + decalage);
                pts.push_back(y2);
            }
        }

        if (a != 0)
        {
            float x1 = -cc / a;
            if (x1 >= 0 && x1 <= w)
            {
                pts.push_back(x1 + decalage);
                pts.push_back(0);
            }

            float x2 = (-cc - b * h) / a;
            if (x2 >= 0 && x2 <= w)
            {
                pts.push_back(x2 + decalage);
                pts.push_back(h);
            }
        }

        if (pts.size() >= 4)
            drawLine(pts[0], pts[1], pts[2], pts[3], c);
    }
}

int main(int argc, char *argv[])
{
    srand((unsigned int)time(0));

    string s1 = argc > 1 ? argv[1] : srcPath("im1.jpg");
    string s2 = argc > 2 ? argv[2] : srcPath("im2.jpg");

    // Load and display images
    Image<Color, 2> I1, I2;
    if (!load(I1, s1) ||
        !load(I2, s2))
    {
        cerr << "Unable to load images" << endl;
        return 1;
    }
    int w = I1.width();
    openWindow(2 * w, I1.height());
    display(I1, 0, 0);
    display(I2, w, 0);

    vector<Match> matches;
    algoSIFT(I1, I2, matches);
    const int n = (int)matches.size();
    cout << " matches: " << n << endl;
    drawString(100, 20, std::to_string(n) + " matches", RED);
    click();

    FMatrix<float, 3, 3> F = computeF(matches);
    cout << "F=" << endl
         << F;

    // Redisplay with matches
    display(I1, 0, 0);
    display(I2, w, 0);
    for (size_t i = 0; i < matches.size(); i++)
    {
        Color c(rand() % 256, rand() % 256, rand() % 256);
        fillCircle(matches[i].x1 + 0, matches[i].y1, 2, c);
        fillCircle(matches[i].x2 + w, matches[i].y2, 2, c);
    }
    drawString(100, 20, to_string(matches.size()) + "/" + to_string(n) + " inliers", RED);
    click();

    // Redisplay without SIFT points
    display(I1, 0, 0);
    display(I2, w, 0);
    displayEpipolar(I1, I2, F);

    endGraphics();
    return 0;
}
