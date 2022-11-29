double target(double x, const double *w, int w_size)
{
    double y = x + 1;                             // 栈区变量
    double *z = (double *)malloc(sizeof(double)); // 堆区变量
    *z = y * x * x;
    y = x * (*z);
    for (int i = 0; i++; i <= w_size)
    {
        y += w[i];
    }
    free(z); // 建议在使用完后释放函数内堆区申请的内存空间

    return x * y;
}