#ifndef __RAY_MATH_TOOLKIT_H
#define __RAY_MATH_TOOLKIT_H

#include <math.h>
#include <stdio.h>
#include <assert.h>

static void normalize(double *v)
{
    const double a = 1.0;

    asm("movsd  (%rdi),%xmm1;"
        "movsd  0x8(%rdi),%xmm2;"
        "movsd  0x10(%rdi),%xmm3;"
        "mulsd  %xmm1,%xmm1;"
        "mulsd  %xmm2,%xmm2;"
        "mulsd  %xmm3,%xmm3;"
        "addsd  %xmm2,%xmm1;"
        "addsd  %xmm3,%xmm1;"
        "sqrtsd %xmm1,%xmm1;"
        "divsd  %xmm1,%xmm0;"
        "movsd  (%rdi),%xmm1;"
        "movsd  0x8(%rdi),%xmm2;"
        "movsd  0x10(%rdi),%xmm3;"
        "mulsd  %xmm0,%xmm1;"
        "mulsd  %xmm0,%xmm2;"
        "mulsd  %xmm0,%xmm3;"
        "movsd  %xmm1, (%rdi);"
        "movsd  %xmm2, 0x8(%rdi);"
        "movsd  %xmm3, 0x10(%rdi);"
       );
}

static double length(const double *v)
{
    asm("movsd  (%rdi),%xmm1;"
        "movsd  0x8(%rdi),%xmm0;"
        "mulsd  %xmm1,%xmm1;"
        "mulsd  %xmm0,%xmm0;"
        "movsd  0x10(%rdi),%xmm2;"
        "mulsd  %xmm2,%xmm2;"
        "addsd  %xmm1,%xmm0;"
        "addsd  %xmm2,%xmm0;"
        "sqrtsd %xmm0,%xmm0;");
}

static void add_vector(const double *a, const double *b, double *out)
{
    asm("movsd  (%rdi),%xmm0;"
        "addsd  (%rsi),%xmm0;"
        "movsd  %xmm0,(%rdx);"
        "movsd  0x8(%rdi),%xmm0;"
        "addsd  0x8(%rsi),%xmm0;"
        "movsd  %xmm0,0x8(%rdx);"
        "movsd  0x10(%rdi),%xmm0;"
        "addsd  0x10(%rsi),%xmm0;"
        "movsd  %xmm0,0x10(%rdx);");
}

static void subtract_vector(const double *a, const double *b, double *out)
{
    asm("movsd  (%rdi),%xmm0;"
        "subsd  (%rsi),%xmm0;"
        "movsd  %xmm0,(%rdx);"
        "movsd  0x8(%rdi),%xmm0;"
        "subsd  0x8(%rsi),%xmm0;"
        "movsd  %xmm0,0x8(%rdx);"
        "movsd  0x10(%rdi),%xmm0;"
        "subsd  0x10(%rsi),%xmm0;"
        "movsd  %xmm0,0x10(%rdx);");
}

static void multiply_vectors(const double *a, const double *b, double *out)
{
    asm("movsd  (%rdi),%xmm0;"
        "mulsd  (%rsi),%xmm0;"
        "movsd  %xmm0,(%rdx);"
        "movsd  0x8(%rdi),%xmm0;"
        "mulsd  0x8(%rsi),%xmm0;"
        "movsd  %xmm0,0x8(%rdx);"
        "movsd  0x10(%rdi),%xmm0;"
        "mulsd  0x10(%rsi),%xmm0;"
        "movsd  %xmm0,0x10(%rdx);");
}

static void multiply_vector(const double *a, double b, double *out)
{
    asm("movsd  (%rdi),%xmm1;"
        "mulsd  %xmm0,%xmm1;"
        "movsd  %xmm1,(%rsi);"
        "movsd  0x8(%rdi),%xmm1;"
        "mulsd  %xmm0,%xmm1;"
        "movsd  %xmm1,0x8(%rsi);"
        "mulsd  0x10(%rdi),%xmm0;"
        "movsd  %xmm0,0x10(%rsi);");
}

static void cross_product(const double *v1, const double *v2, double *out)
{
    asm("movsd  0x8(%rdi),%xmm0;"
        "movsd  0x10(%rdi),%xmm1;"
        "mulsd  0x10(%rsi),%xmm0;"
        "mulsd  0x8(%rsi),%xmm1;"
        "subsd  %xmm1,%xmm0;"
        "movsd  %xmm0,(%rdx);"
        "movsd  0x10(%rdi),%xmm0;"
        "movsd  (%rdi),%xmm1;"
        "mulsd  (%rsi),%xmm0;"
        "mulsd  0x10(%rsi),%xmm1;"
        "subsd  %xmm1,%xmm0;"
        "movsd  %xmm0,0x8(%rdx);"
        "movsd  (%rdi),%xmm0;"
        "movsd  0x8(%rdi),%xmm1;"
        "mulsd  0x8(%rsi),%xmm0;"
        "mulsd  (%rsi),%xmm1;"
        "subsd  %xmm1,%xmm0;"
        "movsd  %xmm0,0x10(%rdx);");
}

static double dot_product(const double *v1, const double *v2)
{
    asm("movsd  (%rdi),%xmm0;"
        "movsd  0x8(%rdi),%xmm1;"
        "mulsd  (%rsi),%xmm0;"
        "mulsd  0x8(%rsi),%xmm1;"
        "addsd  %xmm1,%xmm0;"
        "movsd  0x10(%rdi),%xmm1;"
        "mulsd  0x10(%rsi),%xmm1;"
        "addsd  %xmm1,%xmm0;");
}

static void scalar_triple_product(const double *u, const double *v, const double *w,
                                  double *out)
{
    asm("movsd  0x8(%rsi),%xmm2;"
        "movsd  0x10(%rsi),%xmm0;"
        "mulsd  0x10(%rdx),%xmm2;"
        "mulsd  0x8(%rdx),%xmm0;"
        "subsd  %xmm0,%xmm2;"
        "movsd  %xmm2,(%rcx);"
        "movsd  0x10(%rsi),%xmm1;"
        "movsd  (%rsi),%xmm0;"
        "mulsd  (%rdx),%xmm1;"
        "mulsd  0x10(%rdx),%xmm0;"
        "subsd  %xmm0,%xmm1;"
        "movsd  %xmm1,0x8(%rcx);"
        "movsd  (%rsi),%xmm0;"
        "movsd  0x8(%rsi),%xmm3;"
        "mulsd  0x8(%rdx),%xmm0;"
        "mulsd  (%rdx),%xmm3;"
        "subsd  %xmm3,%xmm0;"
        "movsd  %xmm0,0x10(%rcx);"
        "mulsd  (%rdi),%xmm2;"
        "movsd  %xmm2,(%rcx);"
        "mulsd  0x8(%rdi),%xmm1;"
        "movsd  %xmm1,0x8(%rcx);"
        "mulsd  0x10(%rdi),%xmm0;"
        "movsd  %xmm0,0x10(%rcx);");
}

static double scalar_triple(const double *u, const double *v, const double *w)
{
    asm("movsd  0x8(%rdx),%xmm3;"
        "movsd  0x10(%rdx),%xmm2;"
        "movsd  0x10(%rdi),%xmm4;"
        "movapd %xmm3,%xmm0;"
        "movsd  0x8(%rdi),%xmm1;"
        "movsd  (%rdi),%xmm5;"
        "movapd %xmm2,%xmm7;"
        "movsd  (%rdx),%xmm6;"
        "mulsd  %xmm1,%xmm7;"
        "mulsd  %xmm4,%xmm0;"
        "mulsd  %xmm5,%xmm2;"
        "mulsd  %xmm6,%xmm4;"
        "mulsd  %xmm6,%xmm1;"
        "subsd  %xmm7,%xmm0;"
        "mulsd  %xmm5,%xmm3;"
        "subsd  %xmm4,%xmm2;"
        "mulsd  (%rsi),%xmm0;"
        "subsd  %xmm3,%xmm1;"
        "mulsd  0x8(%rsi),%xmm2;"
        "mulsd  0x10(%rsi),%xmm1;"
        "addsd  %xmm2,%xmm0;"
        "addsd  %xmm1,%xmm0;");
}

#endif
