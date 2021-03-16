#ifndef ELEMENTO_H
#define ELEMENTO_H

#include <string>
#include <vector>
using namespace std;

class Elemento
{
public:
    Elemento(string, int, int);
    virtual ~Elemento() {}
    virtual void addEstampa(vector<vector<double>> &) = 0;
    int na, nb, nc, nd, jx, jy;
    double valor;
    string nome;
};

// ************************ RESISTOR ***************************
class Resistor : public Elemento
{
public:
    Resistor(string, int, int, double);
    virtual void addEstampa(vector<vector<double>> &);
};

// ************************ TRANSCONDUTOR **********************
class Transcondutor : public Elemento
{
public:
    Transcondutor(string, int, int, int, int, double);
    virtual void addEstampa(vector<vector<double>> &);
};

// ************************ AMP DE TENSÃO **********************
class AmpTensao : public Elemento
{
public:
    AmpTensao(string, int, int, int, int, double);
    virtual void addEstampa(vector<vector<double>> &);
};

//  *********************** AMP DE CORRENTE ********************
class AmpCorrente : public Elemento
{
public:
    AmpCorrente(string, int, int, int, int, double);
    virtual void addEstampa(vector<vector<double>> &);
};

// ************************ TRANSRESISTOR ***********************
class Transresistor : public Elemento
{
public:
    Transresistor(string, int, int, int, int, double);
    virtual void addEstampa(vector<vector<double>> &);
};

// ************************ FONTE DE CORRENTE *******************
class FonteCorrente : public Elemento
{
public:
    FonteCorrente(string, int, int, const double &, vector<string>, double &);
    virtual void addEstampa(vector<vector<double>> &);

private:
    void setParam(string);
    void setValor(string tipo);
    double g;
    string tipo;
    vector<string> parametros;
    double &dt;
    double valor_anterior;
    const double &tempo;
    double A0, A, f, ta, al, k, c;
    double k2, A2, ta2, ts, td, t1, p, c2;
};

// ************************ FONTE DE TENSÃO *********************
class FonteTensao : public Elemento
{
public:
    FonteTensao(string, int, int, const double &, vector<string>, double &);
    virtual void addEstampa(vector<vector<double>> &);

private:
    void setParam(string);
    void setValor(string);
    bool inicial;
    string tipo;
    vector<string> parametros;
    double &dt;
    double g;
    double valor_anterior;
    const double &tempo;
    double A0, A, f, ta, al, k, c;
    double k2, A2, ta2, ts, td, t1, p, c2;
};

// ************************ AMP OP ******************************
class AmpOp : public Elemento
{
public:
    AmpOp(string, int, int, int, int);
    virtual void addEstampa(vector<vector<double>> &);
};

class TransformadorIdeal : public Elemento
{
public:
    TransformadorIdeal(string, int, int, int, int, double);
    virtual void addEstampa(vector<vector<double>> &);
};

// *********************** CAPACITOR ****************************
class Capacitor : public Elemento
{
public:
    Capacitor(string, int, int, double, double &, vector<double> &);
    virtual void addEstampa(vector<vector<double>> &);

private:
    vector<double> &vi;
    double ji, k, g_anterior;
    double valor, vAnterior;
    double &dt;
    bool inicial1, inicial2, passo1;
};

// ********************** INDUTOR *******************************
class Indutor : public Elemento
{
public:
    Indutor(string, int, int, double, double &, vector<double> &);
    virtual void addEstampa(vector<vector<double>> &);

private:
    vector<double> &vi;
    double vv, g_anterior, k;
    double valor, jAnterior;
    double &dt;
    bool inicial1, inicial2, passo1;
};

// ********************** RESISTOR LINEAR POR PARTES ************
class ResistorLinearPartes : public Elemento
{
public:
    ResistorLinearPartes(string, int, int, vector<string>, vector<double> &, double &, vector<bool> &);
    virtual void addEstampa(vector<vector<double>> &);

private:
    void setValor();
    double v1, j1, v2, j2, v3, j3, v4, j4;
    double j0, g0, j0_anterior, g0_anterior, g, j, k, gm_anterior;
    double &gm; // Caso seja utilizado o gmin stepping. Começa com zero.
    vector<double> &vnr;
    vector<bool> &conv;
    bool inicial;
};

// ********************** CHAVE *********************************
class Chave : public Elemento
{
public:
    Chave(string, int, int, int, int, double, double, double, vector<double> &, double &, vector<bool> &);
    virtual void addEstampa(vector<vector<double>> &);

private:
    vector<bool> &conv;
    bool inicial;
    void setValor();
    double g, g_anterior, gon, goff, vref, gm_anterior;
    double &gm;
    vector<double> &va;
};

#endif