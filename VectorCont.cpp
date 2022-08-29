#include "MD.h"
#include "MD_VV.h"
#include "MD_Round_VV.h"

int main()
{
    setlocale(LC_ALL, "rus");
    int method_Number = 0, point_Count = 0, space = 0, space_TYPE = 0;
    double circle_Radius = 0.;
    try
    {
        cout << "Выберите форму пространства" << endl;
        cout << "1 - прямоугольное, 2 - круглое: ";
        cin >> space;

        if (space < 1 || space > 3)
            throw BadSelectedMethodException();

        if (space == 1)
        {
            cout << "\nВыберите метод моделирования молекулярной динамики" << endl;
            cout << "1 - классический метод Верле, 2 - скоростной метод Верле: ";
            cin >> method_Number;

            if (method_Number < 1 || method_Number > 2)
                throw BadSelectedMethodException();

            cout << "\nВыберите тип пространства" << endl;
            cout << "1 - пространство с периодическими граничными условиями, 2 - ограниченное пространство: ";
            cin >> space_TYPE;

            if (space_TYPE < 1 && space_TYPE > 2)
                throw BadSelectedMethodException();

            cout << "\nВведите количество частиц: ";
            cin >> point_Count;

            if (point_Count < 1 || point_Count > 5000)
                throw PointCountOutOfRange();

            if (method_Number == 1)
            {
                MD* MolecularDynamics = new MD(point_Count, .001);
                MolecularDynamics->Initialization();
                MolecularDynamics->velocity_Normalization();
                MolecularDynamics->Evolution_СV(50000, space_TYPE);
            }
            else if (method_Number == 2)
            {
                MD_VV* MolecularDynamics = new MD_VV(point_Count, .001);
                MolecularDynamics->Initialization();
                MolecularDynamics->Velocity_Normalization();
                MolecularDynamics->Evolution_VV(20000, space_TYPE, 1);
            }
        }
        else if (space == 2)
        {
            cout << "\nВведите радиус пространства: ";
            cin >> circle_Radius;

            if (circle_Radius < 1 || circle_Radius > 200)
                throw CircleRadiusOutOfRange();

            cout << "\nВведите количество частиц: ";
            cin >> point_Count;

            if (point_Count < 1 || point_Count > 5000)
                throw PointCountOutOfRange();

            MD_Round_VV* MolecularDynamics = new MD_Round_VV(point_Count, circle_Radius, .001);
            MolecularDynamics->Initialization();
            MolecularDynamics->Velocity_Normalization();
            MolecularDynamics->Evolution_VV(10000);
        }
        else if (space == 3)
        {
            cout << "\nВведите радиус пространства: ";
            cin >> circle_Radius;

            if (circle_Radius < 1 || circle_Radius > 100000)
                throw CircleRadiusOutOfRange();

            cout << "\nВведите количество точек: ";
            cin >> point_Count;

            if (point_Count < 1 || point_Count > 9999999)
                throw PointCountOutOfRange();

            MD_Round_VV* MolecularDynamics = new MD_Round_VV(point_Count, circle_Radius);
            MolecularDynamics->Initialization();
        }
    }
    catch (Exception& e)
    {
        e.ShowMessage();
        return 1;
    }
    return 0;
}