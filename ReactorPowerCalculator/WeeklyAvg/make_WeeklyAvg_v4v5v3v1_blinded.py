#!/usr/bin/env python3

import calendar
import datetime


OUTFILE = 'WeeklyAvg_v4v5v3v1_blinded.txt'
P17BFILE = 'WeeklyAvg_P17B_by_Beda.txt'

# NOTE P17B and Post17 overlap in week 296 (day 2076 in P17B -> no data in day
# 2077 -> day 2078 in Post17).
# Inclusive range:
POST17_WEEKMIN = 297
POST17_WEEKMAX = 468            # Post17 last day is 3276

FRAC_POWER = 1.

FRAC_U235 = 0.64
FRAC_U238 = 0.08
FRAC_PU239 = 0.25
FRAC_PU241 = 0.03


def week_range(week):
    day0 = datetime.datetime(2011, 12, 24)
    def stamp(week):
        delta = datetime.timedelta(weeks=week)
        return calendar.timegm((day0 + delta).timetuple())
    return stamp(week), stamp(week + 1)


def main():
    outfile = open(OUTFILE, 'w')

    outfile.write(open(P17BFILE).read())

    for week in range(POST17_WEEKMIN, POST17_WEEKMAX + 1):
        tstart, tend = week_range(week)
        for core in range(1, 6 + 1):
            line = f'{week} {core} {tstart} {tend} {FRAC_POWER:.6f} 0'
            line += f' {FRAC_U235:.6f} {FRAC_U238:.6f}'
            line += f' {FRAC_PU239:.6f} {FRAC_PU241:.6f}'
            outfile.write(line + '\n')


if __name__ == '__main__':
    main()
