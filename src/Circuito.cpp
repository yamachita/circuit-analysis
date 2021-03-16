/*Elementos aceitos e linhas do netlist:
Resistor:                                R<nome> <nó1> <nó2> <Resistência>
Indutor:                                 L<nome> <nó1> <nó2> <Indutância>
Capacitor:                               C<nome> <nó1> <nó2> <Capacitância>
Fonte de tensão controlada a tensão:     E<nome> <nóV+> <nóV-> <nóv+> <nóv-> <Av>
Fonte de corrente controlada a corrente: F<nome> <nóI+> <nóI-> <nói+> <nói-> <Ai>
Fonte de corrente controlada a tensão:   G<nome> <nóI+> <nóI-> <nóv+> <nóv-> <Gm>
Fonte de tensão controlada a corrente:   H<nome> <nóV+> <nóV-> <nói+> <nói-> <Rm>
Fonte de corrente:                       I<nome> <nó+> <nó-> <parâmetros>
Fonte de tensão:                         V<nome> <nó+> <nó-> <parâmetros>
Amplificador operacional ideal:          O<nome> <nó saída+> <nó saída-> <nó entrada+> <nó entrada->
Resistor linear por partes:              N<nome> <nó+> <nó-> <4 pontos vi ji >
Transformador ideal:                     K<nome> <nó a> <nó b> <nó c> <nó d> <n>
Chave:                                   $<nome> <nó a> <nó b> <nó controle c> <nó controle d> <gon> <goff> <vref>
Comentário: *<comentário>
*/

#include "Circuito.h"
#include "Elemento.h"
#include <iostream>
#include <fstream>
#include <cmath>
#include <iomanip>
//#define DEBUG
#define TOLG 1e-15

Circuit::Circuit()
{
    modo_ = 1;
    num_variaveis_ = num_elementos_ = num_nos_ = 0;
    dt_ = 0;
    tempo_total_ = 0;
    tempo_ = 0;
    variaveis_.push_back("0");
    sp = false;
    linear = true;
    gmin_ = 0;
    passos_pontos_tabela = 0;
}

// Decide qual o tipo de analise sera feita.
void Circuit::runSistema(string nome_arquivo)
{

    nome_arquivo_ = nome_arquivo;

    criarNetList(getFile(nome_arquivo));

    switch (modo_)
    {
    case 1:
        runDCLinear();
        break;
    case 2:
        runDCNLinear();
        break;
    case 3:
        runSPLinear();
        break;
    case 4:
        runSPNLinear();
        break;
    }

    salvarTabela(modo_);
}

// MODOS

void Circuit::runDCLinear()
{
    criarSistema();
    if (resolverSistema())
        throw "Sistema singular";
    salvarSolucaoTempo();
}

void Circuit::runDCNLinear()
{
    const double erro_max = 0.01;
    const int max_it = 5;

    double erro = 0, gmin_anterior = 500;
    int cont_it = 1, i = 0;
    bool sol_final = true, g_stepping = false;
    vector<double> sol_anterior; // ULTIMA SOLUCAO QUE CONVERGIU

    double fator[] = {0, sqrt(10), sqrt(sqrt(10)), sqrt(sqrt(sqrt(10))), sqrt(sqrt(sqrt(sqrt(10))))};

    criarSistema();
    if (resolverSistema())
        throw "Sistema singular";

    salvarSolucaoRN();

    erro = getErro(false);

    while (erro > erro_max || sol_final == false)
    {

        if (cont_it > max_it)
        {

            g_stepping = true;
            sol_final = false;
            cont_it = 0;

            if (i > 4)
                throw "Sistema nao converge.";
            if (i == 0)
            {
                gmin_ = 1.1;
                //sol_atual_NR_ = sol_anterior;
            }
            else
            {
                gmin_ = gmin_anterior / fator[i];
                sol_atual_NR_ = sol_anterior;
            }

            ++i;
        }

#ifdef DEBUG
        cout << "Gmin: " << gmin_ << endl;
#endif

        getErro(g_stepping);

        sol_anterior_NR_ = sol_atual_NR_;
        modSistemaNR();

        if (resolverSistema())
            throw "Sistema singular";

        salvarSolucaoRN();
        erro = getErro(g_stepping);
        ++cont_it;

        if (g_stepping == true && erro <= erro_max)
        {

            if (gmin_ > 1e-12)
            {
#ifdef DEBUG
                cout << "CONVERGIU" << endl;
#endif
                sol_anterior = sol_atual_NR_;
                gmin_anterior = gmin_;
                gmin_ = gmin_ / 10;
                cont_it = 0;
                i = 1;
            }
            else
            {
                sol_final = true;
            }
        }
    }

    tabela_.push_back(sol_atual_NR_);
}

void Circuit::runSPLinear()
{

    tempo_ = 0;
    criarSistema();

    if (resolverSistema())
    {
        //sol_atual_tempo_.assign((num_variaveis_ + 1), 0);
        //tabela_.push_back(sol_atual_tempo_);
        throw "Sistema singular";
    }
    else
    {
        salvarSolucaoTempo();
    }

    for (double i = dt_; i <= (tempo_total_ + dt_); i += dt_)
    {
        tempo_ = i;

        modSistemaTempo();
        if (resolverSistema())
            throw "Sistema singular";
        salvarSolucaoTempo();
    }
}

void Circuit::runSPNLinear()
{

    int contgm = 0;

    const double erro_max = 1e-6;
    const int max_it = 6;

    double erro = 0, gmin_anterior = 1.1;
    int cont_it = 1, i = 0;
    bool sol_final = true, g_stepping = false;
    vector<double> sol_anterior; // ULTIMA SOLUCAO QUE CONVERGIU

    double fator[] = {0, sqrt(10), sqrt(sqrt(10)), sqrt(sqrt(sqrt(10))), sqrt(sqrt(sqrt(sqrt(10))))};

    for (double j = 0; j <= (tempo_total_ + dt_); j += dt_)
    {

        tempo_ = j;
        erro = 0;
        gmin_anterior = 1.1;
        cont_it = 1;
        i = 0;
        gmin_ = 0;
        sol_final = true,
        g_stepping = false;
        convergiu_.assign(variaveis_.size(), true);

        if (j != 0)
            modSistemaTempo();
        else
            criarSistema();

        if (resolverSistema())
            throw "Sistema singular";

        salvarSolucaoRN();

        erro = getErro(g_stepping);

#ifdef DEBUG
        cout << "Tempo: " << tempo_ << " Erro: " << erro << " Sol_at: " << sol_atual_NR_[1] << " Sol_ant: " << sol_anterior_NR_[1] << endl;
#endif

        while (erro > erro_max || sol_final == false)
        {

            if (cont_it > max_it)
            {

                g_stepping = true;
                sol_final = false;
                cont_it = 0;

                if (i > 4)
                    throw "Sistema nao converge.";
                if (i == 0)
                {
                    gmin_ = 1.1;
                    //sol_atual_NR_ = sol_anterior;
                }
                else
                {
                    gmin_ = gmin_anterior / fator[i];
                    sol_atual_NR_ = sol_anterior;
                }

                ++i;
            }

#ifdef DEBUG
            cout << "Gmin: " << gmin_ << endl;
#endif

            getErro(g_stepping);

            sol_anterior_NR_ = sol_atual_NR_;
            modSistemaNR();

            if (resolverSistema())
                throw "Sistema singular";

            salvarSolucaoRN();
            erro = getErro(g_stepping);
            ++cont_it;

            if (g_stepping == true && erro <= erro_max)
            {

                if (gmin_ > 1e-12)
                {
#ifdef DEBUG
                    cout << "CONVERGIU" << endl;
#endif
                    sol_anterior = sol_atual_NR_;
                    gmin_anterior = gmin_;
                    gmin_ = gmin_ / 10;
                    cont_it = 0;
                    i = 1;
                }
                else
                {
                    sol_final = true;
                }
            }
        }

#ifdef DEBUG
        cout << "sol_anterior: " << sol_anterior_NR_[1] << endl;
        cout << "sol: " << sol_atual_NR_[1] << endl;
        cout << "--------------------------------------" << endl;
#endif

        sol_anterior = sol_atual_NR_;
        sol_atual_tempo_ = sol_atual_NR_;
        tabela_.push_back(sol_atual_tempo_);
        sol_anterior_NR_ = sol_atual_NR_;
    }
}

// Retorna o erro max encontrado entre o passo atual e o anterior.
// Se o valor da solucao for menor que um retorna o erro absoluto,
// caso contrário retorna o erro relativo.
double Circuit::getErro(bool m)
{
    double max = -1;
    double erro;
    for (int i = 1; i < variaveis_.size(); i++)
    {

        if (std::abs(sol_atual_NR_[i]) > 1)
        {
            erro = std::abs((sol_atual_NR_[i] - sol_anterior_NR_[i]) / (sol_atual_NR_[i]));
        }
        else
        {
            erro = std::abs(sol_atual_NR_[i] - sol_anterior_NR_[i]);
        }

        if (erro > 1e-6 && m == true)
        {
            convergiu_[i] = false;
        }

        if (erro > max)
            max = erro;
    }

    return max;
}

// AQUISIÇÃO DE DADOS

vector<string> Circuit::getFile(string nome_arquivo)
{

    vector<string> resp;
    string s;

    ifstream arquivo(nome_arquivo, ios::in);

    if (!arquivo)
    {
        throw "Arquivo invalido";
    }

    while (getline(arquivo, s))
    {
        resp.push_back(s);
    }

    arquivo.close();

#ifdef DEBUG
    cout << "Lendo netlist: " << endl;
    cout << "Titulo: ";
    for (int i = 0; i < resp.size(); i++)
    {
        cout << resp[i] << endl;
    }
#endif

    return resp;
}

// MONTAGEM NETLIST INTERNO

void Circuit::criarNetList(vector<string> arquivo)
{

    bool tran = false;

    // Começa no 1 porque no netlist de entrada a primeira linha é o número de nós.
    for (int i = 1; i < arquivo.size(); i++)
    {

        vector<string> linha = split(arquivo[i]);

        if (linha[0] == ".TRAN")
        {
            tran = true;
            tempo_total_ = stod(linha[1]);
            dt_ = stod(linha[2]);
            passos_pontos_tabela = stod(linha[4]);
            dt_ = dt_ / passos_pontos_tabela;
            continue;
        }

        arquivo[i][0] = toupper(arquivo[i][0]);
        char tipo = arquivo[i][0];

        if (tipo == '*')
            continue; // Comentário

        addElemento(tipo, linha);
    }

    if (tran == false)
        throw "Comando invalido";

    num_nos_ = variaveis_.size() - 1;

    addCorrente();

    num_elementos_ = netlist1_.size() + netlist2_.size() + netlist3_.size();
    num_variaveis_ = variaveis_.size() - 1;
    sol_atual_NR_.assign(variaveis_.size(), 0);
    sol_anterior_NR_.assign(variaveis_.size(), 0);
    convergiu_.assign(variaveis_.size(), true);

    if (sp == true && linear == true)
    {
        modo_ = 3;
    }
    else if (sp == false && linear == false)
    {
        modo_ = 2;
    }
    else if (sp == true && linear == false)
    {
        modo_ = 4;
    }
}

void Circuit::addElemento(char tipo, const vector<string> &dat)
{

    int a = num(dat[1]);
    int b = num(dat[2]);
    int c, d = -1;

    switch (tipo)
    {

    case 'R':
        netlist1_.push_back(make_unique<Resistor>(dat[0], a, b, stod(dat[3])));
        break;
    case 'G':
        c = num(dat[3]);
        d = num(dat[4]);
        netlist1_.push_back(make_unique<Transcondutor>(dat[0], a, b, c, d, stod(dat[5])));
        break;
    case 'E':
        c = num(dat[3]);
        d = num(dat[4]);
        netlist1_.push_back(make_unique<AmpTensao>(dat[0], a, b, c, d, stod(dat[5])));
        break;
    case 'F':
        c = num(dat[3]);
        d = num(dat[4]);
        netlist1_.push_back(make_unique<AmpCorrente>(dat[0], a, b, c, d, stod(dat[5])));
        break;
    case 'H':
        c = num(dat[3]);
        d = num(dat[4]);
        netlist1_.push_back(make_unique<Transresistor>(dat[0], a, b, c, d, stod(dat[5])));
        break;
    case 'I':
    {
        auto fc = make_unique<FonteCorrente>(dat[0], a, b, tempo_, copy(dat, 3, dat.size() - 1), dt_);
        if (dat[3] == "SIN" || dat[3] == "PULSE")
            sp = true;
        netlist2_.push_back(move(fc));
        break;
    }
    case 'V':
    {
        auto ft = make_unique<FonteTensao>(dat[0], a, b, tempo_, copy(dat, 3, dat.size() - 1), dt_);
        if (dat[3] == "SIN" || dat[3] == "PULSE")
            sp = true;
        netlist2_.push_back(move(ft));
        break;
    }
    case 'O':
        c = num(dat[3]);
        d = num(dat[4]);
        netlist1_.push_back(make_unique<AmpOp>(dat[0], a, b, c, d));
        break;
    case 'C':
        sp = true;
        netlist2_.push_back(make_unique<Capacitor>(dat[0], a, b, stod(dat[3]), dt_, sol_atual_tempo_));
        break;
    case 'L':
        sp = true;
        netlist2_.push_back(make_unique<Indutor>(dat[0], a, b, stod(dat[3]), dt_, sol_atual_tempo_));
        break;
    case 'N':
        linear = false;
        netlist3_.push_back(make_unique<ResistorLinearPartes>(dat[0], a, b, copy(dat, 3, dat.size() - 1), sol_atual_NR_, gmin_, convergiu_));
        break;
    case 'K':
        c = num(dat[3]);
        d = num(dat[4]);
        netlist1_.push_back(make_unique<TransformadorIdeal>(dat[0], a, b, c, d, stod(dat[5])));
        break;
    case '$':
        linear = false;
        c = num(dat[3]);
        d = num(dat[4]);
        netlist3_.push_back(make_unique<Chave>(dat[0], a, b, c, d, stod(dat[5]), stod(dat[6]), stod(dat[7]), sol_atual_NR_, gmin_, convergiu_));
        break;
    default:
        throw "Elemento desconhecido";
    }
}

void Circuit::addCorrente()
{
    char tipo;

    for (int i = 0; i < netlist1_.size(); i++)
    {
        tipo = netlist1_[i]->nome[0];

        if (tipo == 'E' || tipo == 'F' || tipo == 'O' || tipo == 'K')
        {
            variaveis_.push_back("j" + netlist1_[i]->nome);
            netlist1_[i]->jx = (variaveis_.size() - 1);
        }
        else if (tipo == 'H')
        {
            variaveis_.push_back("jx" + netlist1_[i]->nome);
            variaveis_.push_back("jy" + netlist1_[i]->nome);
            netlist1_[i]->jx = (variaveis_.size() - 2);
            netlist1_[i]->jy = (variaveis_.size() - 1);
        }
    }

    for (int i = 0; i < netlist2_.size(); i++)
    {
        tipo = netlist2_[i]->nome[0];

        if (tipo == 'V' || tipo == 'L')
        {
            variaveis_.push_back("j" + netlist2_[i]->nome);
            netlist2_[i]->jx = (variaveis_.size() - 1);
        }
    }
}

int Circuit::num(string nome)
{
    int i = 0;

    while (i < variaveis_.size())
    {
        if (nome == variaveis_[i])
        {
            return i; // nó já existe
        }
        else
        {
            i++;
        }
    }
    variaveis_.push_back(nome);
    return (variaveis_.size() - 1); // novo nó
}

// MONTAGEM DO SISTEMA

void Circuit::criarSistema()
{
    vector<vector<double>> vec(num_variaveis_ + 1, vector<double>(num_variaveis_ + 2, 0));
    sistema_ = vec;

    for (int i = 0; i < netlist1_.size(); i++)
    {
        netlist1_[i]->addEstampa(sistema_);

#ifdef DEBUG
        cout << "Sistema apos a estampa de " << netlist1_[i]->nome << endl;
        printSistemaInterno();
#endif
    }

    for (int i = 0; i < netlist2_.size(); i++)
    {
        netlist2_[i]->addEstampa(sistema_);

#ifdef DEBUG
        cout << "Sistema apos a estampa de " << netlist2_[i]->nome << endl;
        printSistemaInterno();
#endif
    }

    for (int i = 0; i < netlist3_.size(); i++)
    {
        netlist3_[i]->addEstampa(sistema_);

#ifdef DEBUG
        cout << "Sistema apos a estampa de " << netlist3_[i]->nome << endl;
        printSistemaInterno();
#endif
    }
}

void Circuit::modSistemaTempo()
{
    for (int i = 0; i < netlist2_.size(); i++)
    {
        netlist2_[i]->addEstampa(sistema_);

#ifdef DEBUG
        cout << "Sistema apos a estampa de " << netlist2_[i]->nome << endl;
        printSistemaInterno();
#endif
    }
}

void Circuit::modSistemaNR()
{

    for (int i = 0; i < netlist3_.size(); i++)
    {

        netlist3_[i]->addEstampa(sistema_);

#ifdef DEBUG
        cout << "Sistema apos a estampa de " << netlist3_[i]->nome << endl;
        printSistemaInterno();
#endif
    }
}

// RESOLUÇÃO DO SISTEMA

bool Circuit::resolverSistema()
{
    int i, j, l, a;
    double t, p;

    vector<vector<double>> vec;

    vec = sistema_;

    for (i = 1; i <= num_variaveis_; i++)
    {
        t = 0.0;
        a = i;
        for (l = i; l <= num_variaveis_; l++)
        {
            if (std::abs(vec[l][i]) > std::abs(t))
            {
                a = l;
                t = vec[l][i];
            }
        }
        if (i != a)
        {
            for (l = 1; l <= num_variaveis_ + 1; l++)
            {
                p = vec[i][l];
                vec[i][l] = vec[a][l];
                vec[a][l] = p;
            }
        }
        if (std::abs(t) < TOLG)
        {
            return true;
        }
        for (j = num_variaveis_ + 1; j > 0; j--)
        { /* Basta j>i em vez de j>0 */
            vec[i][j] /= t;
            p = vec[i][j];
            if (p != 0) /* Evita operacoes com zero */
                for (l = 1; l <= num_variaveis_; l++)
                {
                    if (l != i)
                        vec[l][j] -= vec[l][i] * p;
                }
        }
    }
    sistema_resolvido_ = vec;
    return false;
}

// ARMAZENAR DADOS

void Circuit::salvarSolucaoTempo()
{
    vector<double> aux;
    aux.push_back(0);
    for (int i = 1; i <= num_variaveis_; i++)
    {
        aux.push_back(sistema_resolvido_[i][num_variaveis_ + 1]);
    }
    tabela_.push_back(aux);
    sol_atual_tempo_ = aux;
}

void Circuit::salvarSolucaoRN()
{
    vector<double> aux;
    aux.push_back(0);
    for (int i = 1; i <= num_variaveis_; i++)
    {
        aux.push_back(sistema_resolvido_[i][num_variaveis_ + 1]);
    }
    sol_atual_NR_ = aux;
}

void Circuit::salvarTabela(int mod)
{

    string nome = nome_arquivo_.replace(nome_arquivo_.find(".net"), nome_arquivo_.length(), ".tab");

    ofstream saida(nome, ios::out);
    saida << "t";
    for (int k = 1; k < variaveis_.size(); k++)
    {
        saida << " " << variaveis_[k];
    }

    saida << '\n';

    if (modo_ == 1 || modo_ == 2)
    {
        for (double i = 0; i <= tempo_total_; i += (dt_ * passos_pontos_tabela))
        {
            saida << i;
            for (int j = 1; j < tabela_[0].size(); j++)
            {
                saida << " " << tabela_[0][j];
            }
            saida << '\n';
        }
    }
    else
    {
        for (double i = 0; i < tabela_.size(); i = (i + passos_pontos_tabela))
        {
            saida << i * dt_;
            for (int j = 1; j < tabela_[0].size(); j++)
            {
                saida << " " << tabela_[i][j];
            }
            saida << '\n';
        }
    }

    saida.close();
}

// AUXILIAR

string Circuit::toUpperCase(string s)
{
    string resp = "";
    for (int i = 0; i < s.size(); i++)
    {
        resp += toupper(s[i]);
    }

    return resp;
}

vector<string> Circuit::split(string s)
{
    vector<string> ret;
    typedef string::size_type string_size;
    string_size i = 0;

    while (i != s.size())
    {
        while (i != s.size() && isspace(s[i]))
        {
            ++i;
        }

        string_size j = i;

        while (j != s.size() && !isspace(s[j]))
        {
            j++;
        }
        if (i != j)
        {
            ret.push_back(s.substr(i, j - i));
            i = j;
        }
    }

    return ret;
}

vector<string> Circuit::copy(const vector<string> &vec, int inicio, int fim)
{
    vector<string> resp;
    for (int i = inicio; i <= fim; i++)
    {
        resp.push_back(vec[i]);
    }
    return resp;
}

// DEBUG

void Circuit::printSistemaInterno()
{
    for (int k = 1; k <= num_variaveis_; k++)
    {
        for (int j = 1; j <= (num_variaveis_ + 1); j++)
        {
            if (sistema_[k][j] != 0)
            {
                cout.width(12);
                cout << left << showpos << sistema_[k][j];
            }
            else
            {
                cout.width(12);
                cout << left << "...";
            }
        }
        cout << endl;
    }
}

void Circuit::printSolucaoAtual(int i)
{
    cout << "Sistema resolvido: " << endl;
    for (int k = 1; k <= num_variaveis_; k++)
    {
        for (int j = 1; j <= num_variaveis_ + 1; j++)
        {
            if (sistema_resolvido_[k][j] != 0)
            {
                cout.width(12);
                cout << left << showpos << sistema_resolvido_[k][j];
            }
            else
            {
                cout.width(12);
                cout << left << "...";
            }
        }
        cout << endl;
    }
    cout << "Solucao: " << endl;

    for (int i = 1; i <= num_variaveis_; i++)
    {
        cout.width(12);
        cout << left << variaveis_[i] << sistema_resolvido_[i][num_variaveis_ + 1] << endl;
    }
}