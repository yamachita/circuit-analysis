#include "Elemento.h"
#include <cmath>
#include <iostream>
#define PI 3.14159265358979323846
using namespace std;
//#define DEBUG

Elemento::Elemento(string n, int a, int b)
{
    nome = n;
    na = a;
    nb = b;
    nc = nd = jx = jy = -1;
    valor = 0;
}

// ************************ RESISTOR ***************************
Resistor::Resistor(string n, int a, int b, double v) : Elemento(n, a, b)
{
    valor = v;
}

void Resistor::addEstampa(vector<vector<double>> &sistema)
{
    double g = 1 / valor;
    sistema[na][na] += g;
    sistema[nb][nb] += g;
    sistema[na][nb] -= g;
    sistema[nb][na] -= g;
}

// ************************ TRANSCONDUTOR **********************
Transcondutor::Transcondutor(string n, int a, int b, int c, int d, double v) : Elemento(n, a, b)
{
    nc = c;
    nd = d;
    valor = v;
}

void Transcondutor::addEstampa(vector<vector<double>> &sistema)
{
    double g = valor;
    sistema[na][nc] += g;
    sistema[nb][nd] += g;
    sistema[na][nd] -= g;
    sistema[nb][nc] -= g;
}

// ************************ AMP DE TENSÃO **********************
AmpTensao::AmpTensao(string n, int a, int b, int c, int d, double v) : Elemento(n, a, b)
{
    nc = c;
    nd = d;
    valor = v;
}

void AmpTensao::addEstampa(vector<vector<double>> &sistema)
{
    double g = valor;
    sistema[na][jx] += 1;
    sistema[nb][jx] -= 1;
    sistema[jx][na] -= 1;
    sistema[jx][nb] += 1;
    sistema[jx][nc] += g;
    sistema[jx][nd] -= g;
}

//  *********************** AMP DE CORRENTE ********************
AmpCorrente::AmpCorrente(string n, int a, int b, int c, int d, double v) : Elemento(n, a, b)
{
    nc = c;
    nd = d;
    valor = v;
}

void AmpCorrente::addEstampa(vector<vector<double>> &sistema)
{
    double g = valor;
    sistema[na][jx] += g;
    sistema[nb][jx] -= g;
    sistema[nc][jx] += 1;
    sistema[nd][jx] -= 1;
    sistema[jx][nc] -= 1;
    sistema[jx][nd] += 1;
}

// ************************ TRANSRESISTOR ***********************
Transresistor::Transresistor(string n, int a, int b, int c, int d, double v) : Elemento(n, a, b)
{
    nc = c;
    nd = d;
    valor = v;
}

void Transresistor::addEstampa(vector<vector<double>> &sistema)
{
    double g = valor;
    sistema[na][jy] += 1;
    sistema[nb][jy] -= 1;
    sistema[nc][jx] += 1;
    sistema[nd][jx] -= 1;
    sistema[jy][na] -= 1;
    sistema[jy][nb] += 1;
    sistema[jx][nc] -= 1;
    sistema[jx][nd] += 1;
    sistema[jy][jx] += g;
}

// *********************** AMP OP *******************************
AmpOp::AmpOp(string n, int a, int b, int c, int d) : Elemento(n, a, b)
{
    nc = c;
    nd = d;
}

void AmpOp::addEstampa(vector<vector<double>> &sistema)
{
    sistema[na][jx] += 1;
    sistema[nb][jx] -= 1;
    sistema[jx][nc] += 1;
    sistema[jx][nd] -= 1;
}

TransformadorIdeal::TransformadorIdeal(string nn, int a, int b, int c, int d, double n) : Elemento(nn, a, b)
{
    nc = c;
    nd = d;
    valor = n;
}

void TransformadorIdeal::addEstampa(vector<vector<double>> &sistema)
{
    sistema[jx][na] += valor;
    sistema[jx][nb] -= valor;
    sistema[jx][nc] -= 1;
    sistema[jx][nd] += 1;
    sistema[na][jx] -= valor;
    sistema[nb][jx] += valor;
    sistema[nc][jx] += 1;
    sistema[nd][jx] -= 1;
}

// ******************************************************************************************************************************

// ************************ FONTE DE CORRENTE *******************
FonteCorrente::FonteCorrente(string n, int a, int b, const double &t, vector<string> param, double &d)
    : Elemento(n, a, b), tempo(t), dt(d)
{
    parametros = param;
    tipo = parametros[0];
    g = 0;
    valor_anterior = 0;
    setParam(tipo);
}

void FonteCorrente::addEstampa(vector<vector<double>> &sistema)
{
    int s = sistema[0].size() - 1;
    setValor(tipo);

    g = valor - valor_anterior;
    valor_anterior = valor;

    sistema[na][s] -= g;
    sistema[nb][s] += g;
}

void FonteCorrente::setValor(string tipo)
{
    if (tipo == "DC")
    {
        valor = stod(parametros[1]);
    }
    else if (tipo == "SIN")
    {
        double tf = (c / f) + ta;
        double t = tempo;

        if (tempo <= ta)
        {
            t = ta;
        }
        else if (tempo > tf)
        {
            t = tf;
        }

        valor = A0 + A * exp(-al * (t - ta)) * sin(2 * PI * f * (t - ta) + PI * k / 180);
    }
    else if (tipo == "PULSE")
    {
        double tf = (c2 * p) + ta2;
        double m, b;
        double t = fmod(((tempo)-ta2), p);

        if (ts == 0)
            ts = dt;
        if (td == 0)
            td = dt;

        if (tempo > ta2 && tempo <= tf)
        {

            if (t >= 0 && t < ts)
            {

                m = (A2 - k2) / ts;
                b = k2;
                valor = m * t + b;
            }
            else if (t >= ts && t <= (ts + t1))
            {

                valor = A2;
            }
            else if (t > (ts + t1) && t <= (ts + t1 + td))
            {
                m = ((k2 - A2) / td);
                b = -1 * (m) * (ts + t1) + A2;
                valor = m * t + b;
            }
            else if (t > (ts + t1 + td))
            {

                valor = k2;
            }
        }
        else
        {
            valor = k2;
        }
    }
}

void FonteCorrente::setParam(string tipo)
{
    if (tipo == "SIN")
    {
        A0 = stod(parametros[1]);
        A = stod(parametros[2]);
        f = stod(parametros[3]);
        ta = stod(parametros[4]);
        al = stod(parametros[5]);
        k = stod(parametros[6]);
        c = stod(parametros[7]);
    }
    else if (tipo == "PULSE")
    {
        k2 = stod(parametros[1]);
        A2 = stod(parametros[2]);
        ta2 = stod(parametros[3]);
        ts = stod(parametros[4]);
        td = stod(parametros[5]);
        t1 = stod(parametros[6]);
        p = stod(parametros[7]);
        c2 = stod(parametros[8]);
    }
}

// ************************ FONTE DE TENSÃO *********************
FonteTensao::FonteTensao(string n, int a, int b, const double &t, vector<string> param, double &d)
    : Elemento(n, a, b), tempo(t), dt(d)
{
    parametros = param;
    tipo = parametros[0];
    g = 0;
    inicial = true;
    valor_anterior = 0;
    setParam(tipo);
}

void FonteTensao::addEstampa(vector<vector<double>> &sistema)
{

    int s = sistema[0].size() - 1;
    setValor(tipo);
    g = (valor - valor_anterior);
    valor_anterior = valor;

    if (inicial == false)
    {
        sistema[jx][s] -= g;
    }
    else
    {
        sistema[na][jx] += 1;
        sistema[nb][jx] -= 1;
        sistema[jx][na] -= 1;
        sistema[jx][nb] += 1;
        sistema[jx][s] -= g;
        inicial = false;
    }
}

void FonteTensao::setValor(string tipo)
{
    if (tipo == "DC")
    {
        valor = stod(parametros[1]);
    }
    else if (tipo == "SIN")
    {
        double tf = (c / f) + ta;
        double t = tempo;

        if (tempo <= ta)
        {
            t = ta;
        }
        else if (tempo > tf)
        {
            t = tf;
        }

        valor = A0 + A * exp(-al * (t - ta)) * sin(2 * PI * f * (t - ta) + PI * k / 180);
    }
    else if (tipo == "PULSE")
    {

        double tf = (c2 * p) + ta2;
        double m, b;
        double t = fmod(((tempo)-ta2), p);

        if (ts == 0)
            ts = dt;
        if (td == 0)
            td = dt;

        if ((tempo) > ta2 && (tempo) <= tf)
        {
            if (t >= 0 && t < ts)
            {
                m = (A2 - k2) / ts;
                b = k2;
                valor = m * t + b;
            }
            else if (t >= ts && t <= (ts + t1))
            {
                valor = A2;
            }
            else if (t > (ts + t1) && t <= (ts + t1 + td))
            {
                m = ((k2 - A2) / td);
                b = -1 * (m) * (ts + t1) + A2;
                valor = m * t + b;
            }
            else if (t > (ts + t1 + td))
            {
                valor = k2;
            }
        }
        else if (tempo <= ta2 || tempo > tf)
        {
            valor = k2;
        }
    }
}

void FonteTensao::setParam(string tipo)
{
    if (tipo == "SIN")
    {
        A0 = stod(parametros[1]);
        A = stod(parametros[2]);
        f = stod(parametros[3]);
        ta = stod(parametros[4]);
        al = stod(parametros[5]);
        k = stod(parametros[6]);
        c = stod(parametros[7]);
    }
    else if (tipo == "PULSE")
    {
        k2 = stod(parametros[1]);
        A2 = stod(parametros[2]);
        ta2 = stod(parametros[3]);
        ts = stod(parametros[4]);
        td = stod(parametros[5]);
        t1 = stod(parametros[6]);
        p = stod(parametros[7]);
        c2 = stod(parametros[8]);
    }
}

// *********************** CAPACITOR ****************************
Capacitor::Capacitor(string n, int a, int b, double v, double &d, vector<double> &var)
    : Elemento(n, a, b), vi(var), dt(d)
{
    valor = v;
    inicial1 = true;
    inicial2 = true;
    passo1 = true;
    k = 0;
    g_anterior = 0;
}

void Capacitor::addEstampa(vector<vector<double>> &sistema)
{
    if (inicial1 == false)
    {
        if (inicial2 == true)
        {
            double t = (1 / pow(10, 9));

            sistema[na][na] += (((2 * valor) / dt) - t);
            sistema[nb][nb] += (((2 * valor) / dt) - t);
            sistema[na][nb] -= (((2 * valor) / dt) - t);
            sistema[nb][na] -= (((2 * valor) / dt) - t);
            inicial2 = false;
        }

        double vAtual = vi[na] - vi[nb];

        if (passo1 == false)
        {
            ji = (((2 * valor) / dt) * (vAtual - vAnterior)) - ji;
        }
        else
        {
            ji = vAtual / pow(10, 9);

            passo1 = false;
        }

        double g = (((2 * valor * vAtual) / dt) + ji);
        k = g - g_anterior;
        g_anterior = g;

        sistema[na][sistema[0].size() - 1] += k;
        sistema[nb][sistema[0].size() - 1] -= k;

        vAnterior = vAtual;
    }
    else
    {
        double t = pow(10, 9);
        sistema[na][na] += 1 / t;
        sistema[nb][nb] += 1 / t;
        sistema[na][nb] -= 1 / t;
        sistema[nb][na] -= 1 / t;

        inicial1 = false;
    }
}

// ********************** INDUTOR *******************************
Indutor::Indutor(string n, int a, int b, double v, double &d, vector<double> &var)
    : Elemento(n, a, b), vi(var), dt(d)
{
    valor = v;
    inicial1 = true;
    inicial2 = true;
    passo1 = true;
    k = 0;
    g_anterior = 0;
}

void Indutor::addEstampa(vector<vector<double>> &sistema)
{

    if (inicial1 == false)
    {
        if (inicial2 == true)
        {
            sistema[jx][jx] += ((2 * valor) / dt) - (1 / pow(10, 9));

            inicial2 = false;
        }

        double jAtual = vi[jx];

        if (passo1 == false)
        {
            vv = ((2 * valor) / dt) * (jAtual - jAnterior) - vv;
        }
        else
        {
            vv = jAtual / pow(10, 9);

            passo1 = false;
        }

        double g = (((2 * valor * jAtual) / dt) + vv);
        k = g - g_anterior;
        g_anterior = g;

        sistema[jx][sistema[0].size() - 1] += k;

        jAnterior = jAtual;
    }
    else
    {
        sistema[na][jx] += 1;
        sistema[nb][jx] -= 1;
        sistema[jx][na] -= 1;
        sistema[jx][nb] += 1;
        sistema[jx][jx] += (1 / pow(10, 9));

        inicial1 = false;
    }
}

// ********************** RESISTOR LINEAR POR PARTES ************

ResistorLinearPartes::ResistorLinearPartes(string n, int a, int b, vector<string> pontos, vector<double> &varnr,
                                           double &gg, vector<bool> &cc) : Elemento(n, a, b), vnr(varnr), gm(gg), conv(cc)
{
    v1 = stod(pontos[0]);
    j1 = stod(pontos[1]);
    v2 = stod(pontos[2]);
    j2 = stod(pontos[3]);
    v3 = stod(pontos[4]);
    j3 = stod(pontos[5]);
    v4 = stod(pontos[6]);
    j4 = stod(pontos[7]);
    j0 = 0;
    g0 = 0;
    g = 0;
    j = 0;
    j0_anterior = 0;
    g0_anterior = 0;
    gm_anterior = 0;
    k = 0;
    inicial = true;
}

void ResistorLinearPartes::addEstampa(vector<vector<double>> &sistema)
{

    setValor();

    g = g0 - g0_anterior;
    j = j0 - j0_anterior;

    if (conv[na] == false || conv[nb] == false)
    {
        k = gm - gm_anterior;
        gm_anterior = gm;
    }
    else
    {
        k = 0;
    }

    g0_anterior = g0;
    j0_anterior = j0;

    sistema[na][na] += (g + k);
    sistema[na][nb] -= (g + k);
    sistema[nb][na] -= (g + k);
    sistema[nb][nb] += (g + k);
    sistema[na][sistema[0].size() - 1] -= j;
    sistema[nb][sistema[0].size() - 1] += j;

#ifdef DEBUG
    cout << "  g0: " << g0 << "  j0: " << j0 << endl;
#endif
}

void ResistorLinearPartes::setValor()
{
    double v;
    if (inicial == false)
    {
        v = vnr[na] - vnr[nb];
    }
    else
    {
        v = 0.1;
        inicial = false;
    }

#ifdef DEBUG
    cout << "v: " << v;
#endif

    if (v < v2)
    {

        g0 = ((j2 - j1) / (v2 - v1));
        j0 = j2 - (g0 * v2);
    }
    else if (v >= v2 && v < v3)
    {

        g0 = ((j3 - j2) / (v3 - v2));
        j0 = j3 - (g0 * v3);
    }
    else if (v >= v3)
    {

        g0 = ((j4 - j3) / (v4 - v3));
        j0 = j4 - (g0 * v4);
    }
}

Chave::Chave(string n, int a, int b, int c, int d, double gn, double gf, double v, vector<double> &vv, double &gg,
             vector<bool> &cc) : Elemento(n, a, b), va(vv), gm(gg), conv(cc)
{
    nc = c;
    nd = d;
    gon = gn;
    goff = gf;
    vref = v;
    g = g_anterior = gm_anterior = 0;
    inicial = true;
}

void Chave::setValor()
{
    double vcd;
    if (inicial == false)
    {
        vcd = va[nc] - va[nd];
    }
    else
    {
        vcd = 0;
        inicial = false;
    }

    if (vcd > vref)
    {
        g = gon;
    }
    else
    {
        g = goff;
    }
}

void Chave::addEstampa(vector<vector<double>> &sistema)
{

    setValor();

    double k = 0;

    double gc = g - g_anterior;

    if (conv[na] == false || conv[nb] == false)
    {
        k = gm - gm_anterior;
        gm_anterior = gm;
    }

    g_anterior = g;

    sistema[na][na] += (gc + k);
    sistema[na][nb] -= (gc + k);
    sistema[nb][na] -= (gc + k);
    sistema[nb][nb] += (gc + k);
}