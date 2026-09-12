"""纯物理容量因子核（仅供 global_bcsd patchify 场站流程使用）。

该模块不包含文件发现、区域筛选或 Slurm 编排；输入和输出均为 numpy 数组，
因此每个 patch 作业可以复用同一套经过验证的风光物理计算。
"""
from __future__ import annotations

import numpy as np
from windpowerlib import WindTurbine

SOLAR_CONSTANT_KW = 1.367
MIN_COS_ZENITH = 1e-4
SYS_COEF = 0.8056
T_STC = 25.0
GAMMA_TEMP = 0.005
C1, C2, C3, C4 = 4.3, 0.943, 0.028, 1.528
TURBINE_TYPE = "GE120/2500"
HUB_HEIGHT = 100.0
REF_HEIGHT = 10.0
POWER_LAW_ALPHA = 1.0 / 7.0
CUT_OUT = 25.0


def _erbs_diffuse_fraction(kt: np.ndarray) -> np.ndarray:
    fd = np.where(kt <= 0.22, 1.0 - 0.09 * kt, np.where(
        kt <= 0.80,
        0.9511 - 0.1604 * kt + 4.388 * kt**2 - 16.638 * kt**3 + 12.336 * kt**4,
        0.165,
    ))
    return np.where(kt <= 0.0, 1.0, np.clip(fd, 0.0, 1.0)).astype(np.float32)


def solar_cos_zenith(doy, hour_decimal_utc, lats, lons) -> np.ndarray:
    doy = np.asarray(doy, dtype=np.float64)
    hour = np.asarray(hour_decimal_utc, dtype=np.float64)
    lats = np.asarray(lats, dtype=np.float64)
    lons = np.asarray(lons, dtype=np.float64)
    b = 2.0 * np.pi * (doy - 1.0) / 365.0
    decl = (0.006918 - 0.399912*np.cos(b) + 0.070257*np.sin(b)
            - 0.006758*np.cos(2*b) + 0.000907*np.sin(2*b)
            - 0.002697*np.cos(3*b) + 0.001480*np.sin(3*b))
    eot = 229.18 * (0.000075 + 0.001868*np.cos(b) - 0.032077*np.sin(b)
                    - 0.014615*np.cos(2*b) - 0.040849*np.sin(2*b))
    angle = np.deg2rad((hour[:, None, None] * 60.0 + eot[:, None, None]
                        + 4.0 * lons[None, None, :]) / 4.0 - 180.0)
    lr = np.deg2rad(lats)
    cos_sza = (np.sin(lr)[None, :, None] * np.sin(decl)[:, None, None]
               + np.cos(lr)[None, :, None] * np.cos(decl)[:, None, None] * np.cos(angle))
    return np.maximum(cos_sza, 0.0).astype(np.float32)


def extraterrestrial_horizontal_irradiance(doy, cos_sza) -> np.ndarray:
    b = 2.0 * np.pi * (np.asarray(doy, dtype=np.float64) - 1.0) / 365.0
    e0 = (1.000110 + 0.034221*np.cos(b) + 0.001280*np.sin(b)
          + 0.000719*np.cos(2*b) + 0.000077*np.sin(2*b))
    return (SOLAR_CONSTANT_KW * e0[:, None, None] * cos_sza).astype(np.float32)


def compute_diffuse_fraction_erbs(ghi_kw, doy, hour_decimal_utc, lats, lons):
    cos_sza = solar_cos_zenith(doy, hour_decimal_utc, lats, lons)
    etr = extraterrestrial_horizontal_irradiance(doy, cos_sza)
    with np.errstate(invalid="ignore", divide="ignore"):
        kt = np.where(etr > 1e-6, ghi_kw / etr, 0.0)
    return _erbs_diffuse_fraction(np.clip(kt, 0.0, 1.0)), cos_sza


def plane_irradiance_two_axis_erbs(ghi_kw, fd, cos_sza):
    dhi = ghi_kw * fd
    beam = np.maximum(ghi_kw - dhi, 0.0)
    with np.errstate(invalid="ignore", divide="ignore"):
        dni = np.where(cos_sza > MIN_COS_ZENITH, beam / cos_sza, 0.0)
    plane = dni + dhi * (1.0 + cos_sza) / 2.0
    plane = np.where(cos_sza > 0.0, plane, 0.0)
    return np.maximum(np.nan_to_num(plane), 0.0).astype(np.float32)


def paper_temperature_coefficient(plane_irradiance_kw, tas_c, wind_speed):
    cell = C1 + C2 * tas_c + C3 * plane_irradiance_kw * 1000.0 - C4 * wind_speed
    return np.maximum(np.nan_to_num(1.0 - GAMMA_TEMP * (cell - T_STC)), 0.0).astype(np.float32)


def compute_solar_cf_chunk(rsds_kw, tas_c, uas, vas, doy, hour_decimal_utc, lats, lons):
    fd, cos_sza = compute_diffuse_fraction_erbs(rsds_kw, doy, hour_decimal_utc, lats, lons)
    plane = plane_irradiance_two_axis_erbs(rsds_kw, fd, cos_sza)
    speed = np.sqrt(np.asarray(uas, np.float32)**2 + np.asarray(vas, np.float32)**2)
    return np.clip(plane * paper_temperature_coefficient(plane, tas_c, speed) * SYS_COEF, 0.0, 1.0).astype(np.float32)


def power_law_ratio(hub_height=HUB_HEIGHT, ref_height=REF_HEIGHT, alpha=POWER_LAW_ALPHA):
    return (float(hub_height) / float(ref_height)) ** float(alpha)


def _get_power_curve_arrays():
    turbine = WindTurbine(turbine_type=TURBINE_TYPE, hub_height=HUB_HEIGHT)
    curve = turbine.power_curve
    ws = np.asarray(curve["wind_speed"], dtype=np.float32)
    power = np.asarray(curve["value"], dtype=np.float32) / 1000.0
    rated = float(turbine.nominal_power) / 1000.0
    order = np.argsort(ws); ws, power = ws[order], power[order]
    keep = ws <= CUT_OUT; ws, power = ws[keep], power[keep]
    if ws.size == 0: raise RuntimeError(f"{TURBINE_TYPE} 功率曲线在切出风速内为空")
    if ws[-1] < CUT_OUT:
        extra = np.arange(ws[-1] + 0.5, CUT_OUT + 0.001, 0.5, dtype=np.float32)
        if extra.size == 0 or extra[-1] < CUT_OUT: extra = np.append(extra, np.float32(CUT_OUT))
        ws = np.concatenate([ws, extra]); power = np.concatenate([power, np.full_like(extra, power[-1])])
    return ws, power, rated


def apply_power_curve(ws_hub, ws_curve, pw_curve, rated_kw):
    speed = np.asarray(ws_hub, dtype=np.float32)
    ws_ext = np.concatenate([[0.0], ws_curve, [CUT_OUT + 0.01]])
    pw_ext = np.concatenate([[0.0], pw_curve, [0.0]])
    out = np.interp(speed, ws_ext, pw_ext)
    out = np.where((speed >= CUT_OUT) | ~np.isfinite(speed), 0.0, out)
    return (out / rated_kw).astype(np.float32)


def compute_wind_cf_chunk(uas, vas, ws_curve, pw_curve, rated_kw, ratio):
    speed = np.sqrt(np.asarray(uas, np.float32)**2 + np.asarray(vas, np.float32)**2) * ratio
    return np.clip(apply_power_curve(speed, ws_curve, pw_curve, rated_kw), 0.0, 1.0)
