/*
expr.c

位置の式 f(x, y, z) (convectionexpr / radiationexpr)。入力の解釈のときに逆ポーランド記法の命令列に
翻訳しておき (構文の誤りはここで入力エラー)、面の中心ごとに評価する。

	式     := 項 (('+' | '-') 項)*
	項     := 単項 (('*' | '/') 単項)*
	単項   := ('-' | '+') 単項 | 累乗
	累乗   := 一次 ('^' 単項)?          (右結合 : 2^-x, 2^3^2 = 2^9)
	一次   := 数 | x | y | z | pi | 関数 '(' 式 (',' 式)? ')' | '(' 式 ')'
	関数   : sin cos tan exp log sqrt abs tanh step (1 引数)、min max pow (2 引数)
	         step(u) = (u >= 0 ? 1 : 0)

評価の演算は書いた順そのもの (左結合) なので、同じ式を awk 等で同じ順に書けば同じ double になる
(検証はそれを使う)。
*/

#include "fem.h"
#include "fem_prototype.h"

enum {
	EX_NUM = 0, EX_X, EX_Y, EX_Z,
	EX_ADD, EX_SUB, EX_MUL, EX_DIV, EX_POW, EX_NEG,
	EX_SIN, EX_COS, EX_TAN, EX_EXP, EX_LOG, EX_SQRT, EX_ABS, EX_TANH, EX_STEP,
	EX_MIN, EX_MAX
};

typedef struct {
	const char *s;
	int pos;
	int err;
	expr_t *e;
} exparse_t;


static void ex_emit(exparse_t *p, int op, double v)
{
	if (p->e->n >= EXPR_MAXOP) {
		p->err = 1;
		return;
	}
	p->e->op[p->e->n] = op;
	p->e->val[p->e->n] = v;
	p->e->n++;
}


static void ex_space(exparse_t *p)
{
	while ((p->s[p->pos] == ' ') || (p->s[p->pos] == '\t')) p->pos++;
}


static void ex_expr(exparse_t *p);
static void ex_unary(exparse_t *p);


static void ex_primary(exparse_t *p)
{
	ex_space(p);
	const char c = p->s[p->pos];
	if (c == '(') {
		p->pos++;
		ex_expr(p);
		ex_space(p);
		if (p->s[p->pos] != ')') { p->err = 1; return; }
		p->pos++;
		return;
	}
	if (((c >= '0') && (c <= '9')) || (c == '.')) {
		char *end = NULL;
		const double v = strtod(&p->s[p->pos], &end);
		if (end == &p->s[p->pos]) { p->err = 1; return; }
		p->pos = (int)(end - p->s);
		ex_emit(p, EX_NUM, v);
		return;
	}
	if (((c >= 'a') && (c <= 'z')) || ((c >= 'A') && (c <= 'Z'))) {
		char name[16];
		int n = 0;
		while ((((p->s[p->pos] >= 'a') && (p->s[p->pos] <= 'z')) || ((p->s[p->pos] >= 'A') && (p->s[p->pos] <= 'Z'))
		      || ((p->s[p->pos] >= '0') && (p->s[p->pos] <= '9'))) && (n < 15)) {
			name[n++] = (char)tolower((int)p->s[p->pos]);
			p->pos++;
		}
		name[n] = '\0';
		if (!strcmp(name, "x")) { ex_emit(p, EX_X, 0); return; }
		if (!strcmp(name, "y")) { ex_emit(p, EX_Y, 0); return; }
		if (!strcmp(name, "z")) { ex_emit(p, EX_Z, 0); return; }
		if (!strcmp(name, "pi")) { ex_emit(p, EX_NUM, PI); return; }
		static const char *f1[] = {"sin", "cos", "tan", "exp", "log", "sqrt", "abs", "tanh", "step"};
		static const int o1[] = {EX_SIN, EX_COS, EX_TAN, EX_EXP, EX_LOG, EX_SQRT, EX_ABS, EX_TANH, EX_STEP};
		static const char *f2[] = {"min", "max", "pow"};
		static const int o2[] = {EX_MIN, EX_MAX, EX_POW};
		int op = -1, nargs = 0;
		for (int k = 0; k < 9; k++) {
			if (!strcmp(name, f1[k])) { op = o1[k]; nargs = 1; }
		}
		for (int k = 0; k < 3; k++) {
			if (!strcmp(name, f2[k])) { op = o2[k]; nargs = 2; }
		}
		if (op < 0) { p->err = 1; return; }
		ex_space(p);
		if (p->s[p->pos] != '(') { p->err = 1; return; }
		p->pos++;
		ex_expr(p);
		if (nargs == 2) {
			ex_space(p);
			if (p->s[p->pos] != ',') { p->err = 1; return; }
			p->pos++;
			ex_expr(p);
		}
		ex_space(p);
		if (p->s[p->pos] != ')') { p->err = 1; return; }
		p->pos++;
		ex_emit(p, op, 0);
		return;
	}
	p->err = 1;
}


static void ex_power(exparse_t *p)
{
	ex_primary(p);
	if (p->err) return;
	ex_space(p);
	if (p->s[p->pos] == '^') {
		p->pos++;
		ex_unary(p);
		ex_emit(p, EX_POW, 0);
	}
}


static void ex_unary(exparse_t *p)
{
	ex_space(p);
	if (p->s[p->pos] == '-') {
		p->pos++;
		ex_unary(p);
		ex_emit(p, EX_NEG, 0);
		return;
	}
	if (p->s[p->pos] == '+') {
		p->pos++;
		ex_unary(p);
		return;
	}
	ex_power(p);
}


static void ex_term(exparse_t *p)
{
	ex_unary(p);
	while (!p->err) {
		ex_space(p);
		const char c = p->s[p->pos];
		if ((c != '*') && (c != '/')) break;
		p->pos++;
		ex_unary(p);
		ex_emit(p, ((c == '*') ? EX_MUL : EX_DIV), 0);
	}
}


static void ex_expr(exparse_t *p)
{
	ex_term(p);
	while (!p->err) {
		ex_space(p);
		const char c = p->s[p->pos];
		if ((c != '+') && (c != '-')) break;
		p->pos++;
		ex_term(p);
		ex_emit(p, ((c == '+') ? EX_ADD : EX_SUB), 0);
	}
}


// 式を翻訳する。戻り値 0 = 正常、1 = 構文の誤り (未知の名前・括弧の対応・余分な文字など)
int expr_compile(const char *s, expr_t *e)
{
	exparse_t p;
	p.s = s;
	p.pos = 0;
	p.err = 0;
	p.e = e;
	e->n = 0;
	ex_expr(&p);
	ex_space(&p);
	if (p.err || (p.s[p.pos] != '\0') || (e->n < 1)) {
		e->n = 0;
		return 1;
	}
	// 評価のスタックの深さ (演算は高々 2 個を 1 個にするので、数と変数の数が上限)
	int depth = 0, dmax = 0;
	for (int i = 0; i < e->n; i++) {
		const int op = e->op[i];
		if (op <= EX_Z) depth++;
		else if ((op >= EX_ADD) && (op <= EX_POW)) depth--;
		else if ((op == EX_MIN) || (op == EX_MAX)) depth--;
		if (depth > dmax) dmax = depth;
	}
	if ((depth != 1) || (dmax > EXPR_MAXSTACK)) {
		e->n = 0;
		return 1;
	}
	return 0;
}


double expr_eval(const expr_t *e, double x, double y, double z)
{
	double st[EXPR_MAXSTACK];
	int sp = 0;
	for (int i = 0; i < e->n; i++) {
		switch (e->op[i]) {
		case EX_NUM: st[sp++] = e->val[i]; break;
		case EX_X: st[sp++] = x; break;
		case EX_Y: st[sp++] = y; break;
		case EX_Z: st[sp++] = z; break;
		case EX_ADD: sp--; st[sp - 1] = st[sp - 1] + st[sp]; break;
		case EX_SUB: sp--; st[sp - 1] = st[sp - 1] - st[sp]; break;
		case EX_MUL: sp--; st[sp - 1] = st[sp - 1] * st[sp]; break;
		case EX_DIV: sp--; st[sp - 1] = st[sp - 1] / st[sp]; break;
		case EX_POW: sp--; st[sp - 1] = pow(st[sp - 1], st[sp]); break;
		case EX_MIN: sp--; st[sp - 1] = ((st[sp] < st[sp - 1]) ? st[sp] : st[sp - 1]); break;
		case EX_MAX: sp--; st[sp - 1] = ((st[sp] > st[sp - 1]) ? st[sp] : st[sp - 1]); break;
		case EX_NEG: st[sp - 1] = -st[sp - 1]; break;
		case EX_SIN: st[sp - 1] = sin(st[sp - 1]); break;
		case EX_COS: st[sp - 1] = cos(st[sp - 1]); break;
		case EX_TAN: st[sp - 1] = tan(st[sp - 1]); break;
		case EX_EXP: st[sp - 1] = exp(st[sp - 1]); break;
		case EX_LOG: st[sp - 1] = log(st[sp - 1]); break;
		case EX_SQRT: st[sp - 1] = sqrt(st[sp - 1]); break;
		case EX_ABS: st[sp - 1] = fabs(st[sp - 1]); break;
		case EX_TANH: st[sp - 1] = tanh(st[sp - 1]); break;
		case EX_STEP: st[sp - 1] = ((st[sp - 1] >= 0) ? 1.0 : 0.0); break;
		default: break;
		}
	}
	return ((sp == 1) ? st[0] : 0);
}
