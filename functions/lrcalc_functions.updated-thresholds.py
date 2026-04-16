from math import exp, log, sqrt
from scipy.integrate import quad
from scipy.stats import norm
import sys


def adaptive_gaussian_quadrature(f, a, b, tol=1e-8):
    result, error = quad(f, a, b, epsabs=tol)
    return result


def binomial_likelihood(n, k, p):
    return (p ** k) * ((1 - p) ** (n - k))


def calculate_odds_ratio(case_var_count, control_var_count, case_nonvar_count, control_nonvar_count, confidence=95):
    haldane = False
    if case_var_count == 0 or control_var_count == 0:
        case_var_count += 0.5
        control_var_count += 0.5
        case_nonvar_count += 0.5
        control_nonvar_count += 0.5
        haldane = True

    odds_ratio = (case_var_count/case_nonvar_count)/(control_var_count/control_nonvar_count)
    sem = sqrt(sum([1/x for x in (case_var_count, control_var_count, case_nonvar_count, control_nonvar_count)]))

    # Calculate coefficient for CI calculation using norm PPF function
    num_sems = norm.ppf(confidence/100 + ((1-confidence/100)/2))
    lower_ci = exp(log(odds_ratio) - num_sems * sem)
    upper_ci = exp(log(odds_ratio) + num_sems * sem)

    return lower_ci, upper_ci


def convert_or_hypotheses_to_p(h0_upper, h1_lower, case_var_count, control_var_count,
                               case_nonvar_count, control_nonvar_count):

    # For H0: if lower bound is None, p0_lower should be 0
    p_h0l = 0

    # Historic implementation of p_h0u
    # p_h0u = (h0_upper * case_nonvar_count) / (h0_upper * case_nonvar_count + control_nonvar_count)

    # Updated implementation using quadratic formula to calculate expected prob. under target OR(non-assoc) = h0_upper
    p_h0u = calc_p_via_quadratic(h0_upper, case_var_count, control_var_count, case_nonvar_count, control_nonvar_count)

    # Historic implementation of p_h1l
    # p_h1l = (h1_lower * case_nonvar_count) / (h1_lower * case_nonvar_count + control_nonvar_count)

    # Updated implementation using quadratic formula to calculate expected prob. under target OR(assoc) = h1_lowe
    p_h1l = calc_p_via_quadratic(h1_lower, case_var_count, control_var_count, case_nonvar_count, control_nonvar_count)

    # For H1: if upper bound is None, p1_upper should be 1
    p_h1u = 1

    p_h0_range = (p_h0l, p_h0u)
    p_h1_range = (p_h1l, p_h1u)

    return p_h0_range, p_h1_range


def calc_p_via_quadratic(target, case_var_count, control_var_count, case_nonvar_count, control_nonvar_count):

    """ Calculating a scaled value of p when using N1/(N1+N2) as base probability requires use of the quadratic
    formula, which is implemented here. """

    if target == 1:
        return (case_var_count + case_nonvar_count) / (
                    case_var_count + case_nonvar_count + control_var_count + control_nonvar_count)

    p2_coeff = (target - 1) * (case_var_count + control_var_count)
    p_coeff = -(target * (
            2 * case_var_count + control_var_count + case_nonvar_count) + control_nonvar_count - case_var_count)
    const_coeff = target * (case_var_count + case_nonvar_count)
    discriminant = sqrt(p_coeff ** 2 - 4 * p2_coeff * const_coeff)

    return (-p_coeff - discriminant) / (2 * p2_coeff)


def calc_lr(case_var_count, total_case_count, control_var_count, total_control_count, PATHOGENIC_OR_THRESHOLD=5, BENIGN_OR_THRESHOLD=1, CONFIDENCE_LEVEL=0):
    #print('Calcing LR')
    case_nonvar_count = (total_case_count) - (case_var_count)
    control_nonvar_count = (total_control_count) - (control_var_count)

    n = (case_var_count) + (control_var_count)
    k = (case_var_count)

    h0_range, h1_range = convert_or_hypotheses_to_p(BENIGN_OR_THRESHOLD, PATHOGENIC_OR_THRESHOLD,
                                                    case_var_count, control_var_count,
                                                    case_nonvar_count, control_nonvar_count)

    binomial_likelihood_fixed = lambda p: binomial_likelihood(n, k, p)

    h0_l = adaptive_gaussian_quadrature(binomial_likelihood_fixed, h0_range[0], h0_range[1])
    h1_l = adaptive_gaussian_quadrature(binomial_likelihood_fixed, h1_range[0], h1_range[1])

    # odds_base = case_nonvar_count / control_nonvar_count
    # odds_benign = BENIGN_OR_THRESHOLD * odds_base
    # odds_path = PATHOGENIC_OR_THRESHOLD * odds_base
    #
    # carrier_count = int(case_var_count) + int(control_var_count)
    #
    # benign_pseudo_cases = odds_benign * carrier_count / (odds_benign + 1)
    # benign_pseudo_controls = carrier_count - benign_pseudo_cases
    #
    # _, benign_uci = \
    #     calculate_odds_ratio(benign_pseudo_cases, benign_pseudo_controls,
    #                          case_nonvar_count, control_nonvar_count,
    #                          confidence=CONFIDENCE_LEVEL)
    #
    # path_pseudo_cases = odds_path * carrier_count / (odds_path + 1)
    # path_pseudo_controls = carrier_count - path_pseudo_cases
    #
    # path_lci, _ = \
    #     calculate_odds_ratio(path_pseudo_cases, path_pseudo_controls,
    #                          case_nonvar_count, control_nonvar_count,
    #                          confidence=CONFIDENCE_LEVEL)
    #
    # benign_uci = float(benign_uci) * odds_base / (float(benign_uci) * odds_base + 1)
    # path_lci = float(path_lci) * odds_base / (float(path_lci) * odds_base + 1)
    #
    # binomial_likelihood_fixed = lambda p: binomial_likelihood(n, k, p)
    #
    # table_h0_l = adaptive_gaussian_quadrature(binomial_likelihood_fixed, 0, benign_uci)
    # table_h1_l = adaptive_gaussian_quadrature(binomial_likelihood_fixed, path_lci, 1)
    #
    # table_total_l = adaptive_gaussian_quadrature(binomial_likelihood_fixed, 0, 1)

    # table_hl_percent = (table_h0_l + table_h1_l) * 100 / table_total_l
    #
    # table_lr = table_h1_l / table_h0_l
    table_lr = h1_l / h0_l

    try:
        acmg_eps = log(table_lr, 2.08)
    except ValueError:
        return table_lr, "NA"

    return table_lr, acmg_eps

# lr, acmg_eps = calc_lr(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4])
# print(f'LR: {lr}')
# print(f'ACMG EPs: {acmg_eps}')
