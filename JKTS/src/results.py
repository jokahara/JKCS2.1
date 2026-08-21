import os
import math
import time
from collections import namedtuple

# Rich result of a single-channel MC-TST rate-constant calculation.
RateResult = namedtuple('RateResult',
    ['k', 'kappa', 'Ea', 'Q_TS', 'Q_reactant', 'sigma', 'n_ts', 'n_reactant', 'imag', 'T'])
RateResult.__new__.__defaults__ = (float('nan'), float('nan'), None, None, None, 1, 0, 0, None, 298.15)


K_UNITS = 'cm3 molecule-1 s-1'

# Column layout of a .rates.tsv row (after the leading channel name). New fields are
# appended so that .rates.tsv files written by older versions still read back.
(C_METHOD, C_SIGMA, C_EA, C_KAPPA, C_K, C_QTS, C_T, C_NOTE,
 C_NTS, C_NREAC, C_IMAG, C_BASIS, C_REACTION, C_QREAC, C_STAMP) = range(15)
N_COLS = 15


def format_rate(result):
    if result.k is None or (isinstance(result.k, float) and math.isnan(result.k)):
        return f"k({result.T} K) could not be computed (missing energies)"
    return (f"k({result.T} K) = {result.k:.3e} {K_UNITS}  "
            f"(kappa = {result.kappa:.2f}, Ea = {result.Ea:.2f} kcal/mol, sigma = {result.sigma})")


def _num(x, fmt):
    if x is None or (isinstance(x, float) and math.isnan(x)):
        return '-'
    return format(x, fmt)


def _get(vals, i):
    return vals[i] if i < len(vals) else ''


def _formula(channel_name):
    return channel_name.split('_H')[0]


def _channel_tag(channel_name):
    i = channel_name.rfind('_H')
    return channel_name[i + 1:] if i != -1 else channel_name


def record_rate(rates_dir, channel_name, result, method=None, note='', basis_set=None, reaction=None):
    store = os.path.join(rates_dir, '.rates.tsv')
    rows = {}
    if os.path.exists(store):
        with open(store) as f:
            for line in f:
                parts = line.rstrip('\n').split('\t')
                if len(parts) >= 7:
                    rows[parts[0]] = parts[1:]
    row = [''] * N_COLS
    row[C_METHOD] = method or ''
    row[C_SIGMA] = str(result.sigma)
    row[C_EA] = _num(result.Ea, '.4f')
    row[C_KAPPA] = _num(result.kappa, '.4f')
    row[C_K] = repr(result.k)
    row[C_QTS] = _num(result.Q_TS, '.6e')
    row[C_T] = str(result.T)
    row[C_NOTE] = note
    row[C_NTS] = str(result.n_ts or '')
    row[C_NREAC] = str(result.n_reactant or '')
    row[C_IMAG] = _num(-abs(result.imag) if result.imag is not None else None, '.1f')
    row[C_BASIS] = basis_set or ''
    row[C_REACTION] = reaction or ''
    row[C_QREAC] = _num(result.Q_reactant, '.6e')
    row[C_STAMP] = time.strftime('%Y-%m-%d %H:%M')
    rows[channel_name] = row
    with open(store, 'w') as f:
        for ch, vals in rows.items():
            f.write('\t'.join([ch] + list(vals)) + '\n')
    _write_summary(os.path.join(rates_dir, 'Rate_constants.txt'), rows)


def _write_summary(path, rows):
    def hnum(ch):
        tag = _channel_tag(ch)
        return int(tag[1:]) if tag[1:].isdigit() else 0

    def kval(ch):
        try:
            return float(_get(rows[ch], C_K))
        except ValueError:
            return float('nan')

    channels = sorted(rows, key=hnum)
    if not channels:
        return

    def first(col):
        # Legacy .rates.tsv rows lack the newer columns, so take the first channel that has one
        return next((_get(rows[ch], col) for ch in channels if _get(rows[ch], col)), '')

    formula = _formula(channels[0])
    level = ' '.join(x for x in (first(C_METHOD), first(C_BASIS)) if x)
    T = first(C_T)
    radical = first(C_REACTION) or 'OH'
    total = sum(k for k in (kval(ch) for ch in channels) if not math.isnan(k))

    W = 98
    head = f' {formula} + {radical}' + (f'   ({level})' if level else '')
    right = f'T = {T} K'
    n_reac = first(C_NREAC)
    stamp = first(C_STAMP)
    sub_left = f' Reactant conformers: {n_reac}' if n_reac else ''
    sub_right = f'written {stamp}' if stamp else ''

    out = ['=' * W, f'{head:<{W - len(right) - 1}}{right}']
    if sub_left or sub_right:
        out.append(f'{sub_left:<{W - len(sub_right) - 1}}{sub_right}')
    out += ['=' * W,
            f" {'Channel':<9}{'sigma':>6}{'N_TS':>6}{'Ea/kcal/mol':>13}{'kappa':>8}"
            f"{'nu_imag/cm-1':>14}{'k / ' + K_UNITS:>32}{'%':>8}",
            ' ' + '-' * (W - 2)]

    notes = {}
    for ch in channels:
        vals = rows[ch]
        note = _get(vals, C_NOTE)
        k = kval(ch)
        kf = 'failed' if math.isnan(k) else f'{k:.3e}'
        branch = '-' if (math.isnan(k) or total <= 0) else f'{100.0 * k / total:.1f}'
        tag = _channel_tag(ch)
        if note:
            notes[tag] = note
            tag += '*'
        out.append(f" {tag:<9}{_get(vals, C_SIGMA):>6}{_get(vals, C_NTS) or '-':>6}"
                   f"{_get(vals, C_EA):>13}{_get(vals, C_KAPPA):>8}{_get(vals, C_IMAG) or '-':>14}"
                   f"{kf:>32}{branch:>8}")

    out += [' ' + '-' * (W - 2),
            f' TOTAL  k({T} K) = {total:.3e} {K_UNITS}', '=' * W]
    for tag, note in notes.items():
        out.append(f' * {tag}: {note}')
    with open(path, 'w') as f:
        f.write('\n'.join(out) + '\n')


def write_molecule_summary(path, molecules, title=None):
    W = 96
    out = ['=' * W, ' Molecule summary' + (f'  ({title})' if title else ''), '=' * W,
           f" {'Name':<34}{'Step':<12}{'E_elec / Ha':>16}{'ZPE / Ha':>14}{'Q':>16}",
           ' ' + '-' * (W - 2)]
    for m in molecules:
        out.append(f" {(m.name or '')[:33]:<34}{str(getattr(m, 'current_step', '') or '')[:11]:<12}"
                   f"{_num(getattr(m, 'electronic_energy', None), '.6f'):>16}"
                   f"{_num(getattr(m, 'zero_point', None), '.6f'):>14}"
                   f"{_num(getattr(m, 'Q', None), '.3e'):>16}")
    out.append('=' * W)
    with open(path, 'w') as f:
        f.write('\n'.join(out) + '\n')
