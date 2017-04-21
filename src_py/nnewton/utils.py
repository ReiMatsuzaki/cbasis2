def cond(q, true_case, false_case):
    if(q):
        return true_case
    else:
        return false_case


def flatten(xs):
    return reduce(lambda a, b: a+b, xs)


def uniq(xs):
    return list(set(xs))


def with_index(xs):
    return zip(xs, range(len(xs)))


def repeat(x, num):
    return [x for i in range(num)]


def replace_at(xs, i, y):
    ys = list(xs)
    ys[i] = y
    return ys


def list_indexed(ix_list, num):
    xs = repeat(0, num)
    for (i, x) in ix_list:
        xs[i] = x
    return xs

