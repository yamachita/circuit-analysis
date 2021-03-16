#include "Circuito.h"
#include <iostream>
int main()
{

    Circuit c;
    cout << "Diga o nome do arquivo: " << endl;
    string nome_arquivo = "";
    getline(cin, nome_arquivo);

    try
    {

        c.runSistema(nome_arquivo);
        cout << "Terminado com sucesso." << endl;
    }
    catch (const char *msg)
    {
        cout << msg << endl;
    }

    cin.get();
}