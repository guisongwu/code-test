// 单位正方形 [0,1]×[0,1] 的四边形网格生成脚本
// 保存为 square_quad.geo 后用 GMSH 打开

// 定义几何参数
L = 1.0;          // 正方形边长
N = 3;           // 每条边的网格分段数（可调整）
alpha = 1/1.414213562;
// 3 5 9

// 定义四个角点
Point(1) = {0, 0, 0};  // 左下角 (x,y,z)
Point(2) = {L, 0, L/alpha};  // 右下角
Point(3) = {L, L, L/alpha};  // 右上角
Point(4) = {0, L, 0};  // 左上角

// 定义四条边
Line(1) = {1, 2};  // 底边
Line(2) = {2, 3};  // 右边
Line(3) = {3, 4};  // 顶边
Line(4) = {4, 1};  // 左边

// 定义线循环构成面
Curve Loop(1) = {1, 2, 3, 4};
Plane Surface(1) = {1};

// 设置网格划分参数
Transfinite Curve {1,3} = N;  // 水平边分段数
Transfinite Curve {2,4} = N;  // 垂直边分段数
Transfinite Surface {1};      // 应用结构化网格

// 确保生成四边形网格
Recombine Surface {1};

// 可选：添加物理组（便于后续FEM分析）
Physical Curve("Bottom") = {1};
Physical Curve("Right") = {2};
Physical Curve("Top") = {3};
Physical Curve("Left") = {4};
Physical Surface("Domain") = {1};

// 设置网格算法为结构化（重要！）
Mesh.RecombinationAlgorithm = 1;  // 1=Blossom, 2=Simple
Mesh.RecombineAll = 1;
Mesh.SubdivisionAlgorithm = 1;    // 1=All Quad
Mesh.Algorithm = 8;               // 8=Frontal/Delaunay for Quad

// 生成2D网格
Mesh 2;

// 保存网格文件（取消注释使用）
// Save "square_quad.msh";
