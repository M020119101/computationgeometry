#include "GeomBase.h"

int main()
{
    //点、向量、矩阵的用法
    Point3D p1(1, 1, 1);
    Point3D p2(2, 2, 2);
    double dis = p1.distance(p2);
    std::cout << p1.getX() - p2.getX() << " " << p1.getY() - p2.getY() << " " << p1.getZ() - p2.getZ() << std::endl;
    Vector3D v = p1 - p2;
    std::cout << v.getDx() << " " << v.getDy() << " " << v.getDz() << std::endl;
    p1.translate(v);
    Point3D p = p2.translated(v);
    Point3D p3 = p2 + v;
    std::cout << p1.getX() << " " << p1.getY() << " " << p1.getZ() << std::endl;
    std::cout << p2.getX() << " " << p2.getY() << " " << p2.getZ() << std::endl;
    std::cout << p.getX() << " " << p.getY() << " " << p.getZ() << std::endl;
    std::cout << p3.getX() << " " << p3.getY() << " " << p3.getZ() << std::endl;
    std::cout << dis << std::endl;

    Matrix3D m = Matrix3D::createTranslateMatrix(1, 2, 3);
    Point3D p4 = p1 * m;
    std::cout << p4.getX() << " " << p4.getY() << " " << p4.getZ() << std::endl;
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            std::cout << m.a[i][j] << " ";
        }
        std::cout << std::endl;
    }

}

