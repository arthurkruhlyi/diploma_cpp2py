# 	int
# 		Indicator(-1),                              // Переменная-условие продолжения расчета.
# 		SeparationSurface(0);                       // Переменная хранящая узел "разделяющий" слои материалов.
# 	const int
# 		K = 1000,                                    // Количество узлов по ширине.
# 		N = 500;                                    // Количество узлов по глубине.
# 		*zn = new double[N + 1], *rk = new double[K + 1]
# 	const double
# 		BeginTemperature = 300,                     // Начальная температура образца (К).
# //////////////////////////////////////////////////////
# 		H = 5e-4,                                   // Глубина образца (м).
# 		R = 11e-4,                                  // Ширина образца (м).
# 		ThicknessOfFirstLayer = 2e-4,               // Толщина внешнего слоя (м).


def mesh_creation(width_nods_amount, height_nodes_amount):
    if width_nods_amount < 0 or height_nodes_amount < 0:
        raise ValueError("Input values can't be negative!")

    if [type(width_nods_amount), type(height_nodes_amount)] != [int, int]:
        raise TypeError("Input values must be int type only")

    mesh = [widths_node for widths_node in range(width_nods_amount)]
    mesh = [[height_nodes for height_nodes in range(height_nodes_amount)] for node in mesh]

    return mesh


if __name__ == '__main__':
    print(mesh_creation(0, 0))

# MeshCreator(H, R, ThicknessOfFirstLayer, K, N, zn, rk, SeparationSurface, Indicator);
#
# // Функция создающая расчетную область.
# void MeshCreator(
#     double const &rH
#   , double const &rR
#   , double const &rThicknessOfFirstLayer
#   , int const &rK
#   , int const &rN
#   , double *zn
#   , double *rk
#   , int &rSeparationSurface
#   , int &rIndicator
#)
# {
# 	int
# 		i, switcher(0);
# 	double
# 		qr, qz, hz, lr, delta, deltaPrev;
# 	const int
# 		qrK_1 = 2,                                  // Параметр прогрессии по ширине.
# 		qzN_1 = 20;                                 // Параметр прогрессии по глубине.
# 	if (rThicknessOfFirstLayer > rH)
# 		rIndicator = 1;
# 	qr = 1.0 * pow(qrK_1, 1.0 / (rK - 1));
# 	qz = 1.0 * pow(qzN_1, 1.0 / (rN - 1));
# 	hz = rH * (1 - qz) / (1 - pow(qz, rN));
# 	lr = rR * (1 - qr) / (1 - pow(qr, rK));
# 	zn[0] = 0;
# 	for (int i(1); i <= rN; i++)
# 		zn[i] = zn[i - 1] + hz * pow(qz, i - 1);
# 	rk[0] = 0;
# 	for (int i(1); i <= rK; i++)
# 		rk[i] = rk[i - 1] + lr * pow(qr, i - 1);
# 	if ((rThicknessOfFirstLayer != 0))                 // Блок который находит границу "раздела материалов".
# 	{
# 		deltaPrev = rThicknessOfFirstLayer;
# 		i = 1;
# 		switcher = 0;
# 		do
# 		{
# 			if (switcher == 0)
# 			{
# 				delta = rThicknessOfFirstLayer - zn[i];
# 				if (deltaPrev > abs(delta))
# 				{
# 					deltaPrev = delta;
# 					i++;
# 				}
# 				else
# 				{
# 					rSeparationSurface = i - 1;
# 					zn[i - 1] = rThicknessOfFirstLayer;
# 					switcher = 1;
# 				}
# 			}
# 			else
# 				break;
# 		} while (i < rK - 1);
# 	}
# }
