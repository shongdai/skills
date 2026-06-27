#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
test_geo_scripts.py — dsh-geo 技能纯函数单元测试
================================================
测试范围：
  - GEO_expMatrix_check.py : classify_folder / EXP_RE / check_single_matrix
  - GEO_smart_sniff.py     : 常量验证（SPECIES_MAP / 平台集合 / NCBI_GENOME_TAG）

运行：
  cd dsh-geo
  python -m pytest tests/test_geo_scripts.py -v
  # 或
  python tests/test_geo_scripts.py
"""

import json
import os
import re
import sys
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

# 将 scripts/ 加入 sys.path 以便导入
_SCRIPTS_DIR = os.path.join(os.path.dirname(__file__), "..", "scripts")
sys.path.insert(0, os.path.abspath(_SCRIPTS_DIR))


# ================================================================
# 从 GEO_expMatrix_check.py 导入纯函数
# ================================================================
from GEO_expMatrix_check import classify_folder, EXP_RE, NCBI_TYPES


class TestClassifyFolder(unittest.TestCase):
    """测试文件夹类型自动识别（6 种类型）。"""

    def test_ncbicount(self):
        """expMatrix 类型含 Count/FPKM/TPM → ncbicount"""
        files = {"expMatrix_Count.csv", "expMatrix_FPKM.csv", "clinical.csv"}
        self.assertEqual(classify_folder(set(files), {"Count", "FPKM"}), "ncbicount")
        self.assertEqual(classify_folder(set(files), {"TPM"}), "ncbicount")

    def test_cel(self):
        """有 .CEL 文件 + probe2gene → cel"""
        files = {"GSMxxx.CEL.gz", "probe2gene.csv", "clinical.csv"}
        self.assertEqual(classify_folder(set(files), set()), "cel")

        # _RAW.tar 也算 cel
        files2 = {"GSE123_RAW.tar", "probe2gene.csv"}
        self.assertEqual(classify_folder(set(files2), set()), "cel")

    def test_probe(self):
        """有 series_matrix + probe2gene → probe"""
        files = {"GSE123_series_matrix.txt.gz", "probe2gene.csv", "clinical.csv"}
        self.assertEqual(classify_folder(set(files), set()), "probe")

    def test_count(self):
        """无 series_matrix 也无 probe2gene，但 exp_types 非空 → count"""
        files = {"expMatrix.csv", "clinical.csv"}
        self.assertEqual(classify_folder(set(files), {"Single"}), "count")

    def test_raw(self):
        """有 series_matrix 但无 expMatrix → raw（已下载未处理）"""
        files = {"GSE123_series_matrix.txt.gz", "clinical.csv"}
        self.assertEqual(classify_folder(set(files), set()), "raw")

    def test_unknown(self):
        """空文件夹 → unknown"""
        self.assertEqual(classify_folder(set(), set()), "unknown")

    def test_ncbi_takes_priority(self):
        """ncbicount 优先级最高：即使有 series_matrix + probe2gene"""
        files = {
            "GSE123_series_matrix.txt.gz",
            "probe2gene.csv",
            "expMatrix_Count.csv",
        }
        self.assertEqual(classify_folder(set(files), {"Count"}), "ncbicount")


class TestExpRegex(unittest.TestCase):
    """测试 expMatrix 文件名正则匹配。"""

    def test_single_platform(self):
        m = EXP_RE.match("expMatrix.csv")
        self.assertIsNotNone(m)
        self.assertEqual(m.group(1), None)

    def test_with_suffix(self):
        m = EXP_RE.match("expMatrix_GPL570.csv")
        self.assertIsNotNone(m)
        self.assertEqual(m.group(1), "GPL570")

    def test_count_type(self):
        m = EXP_RE.match("expMatrix_Count.csv")
        self.assertIsNotNone(m)
        self.assertEqual(m.group(1), "Count")

    def test_no_match(self):
        self.assertIsNone(EXP_RE.match("clinical.csv"))
        self.assertIsNone(EXP_RE.match("probe2gene.csv"))
        self.assertIsNone(EXP_RE.match("expMatrix"))


class TestCheckSingleMatrix(unittest.TestCase):
    """测试单个 expMatrix CSV 的数据质量检查。"""

    @classmethod
    def setUpClass(cls):
        """创建临时 CSV 文件用于测试。"""
        cls.tmpdir = tempfile.mkdtemp(prefix="test_geo_")
        cls.csv_clean = os.path.join(cls.tmpdir, "expMatrix.csv")
        cls.csv_with_na = os.path.join(cls.tmpdir, "expMatrix_WithNA.csv")
        cls.csv_corrupt = os.path.join(cls.tmpdir, "expMatrix_Bad.csv")

        with open(cls.csv_clean, "w") as f:
            f.write("symbol,sample1,sample2\n")
            f.write("TP53,120.5,130.2\n")
            f.write("BRCA1,45.0,50.1\n")

        with open(cls.csv_with_na, "w") as f:
            f.write("symbol,sample1,sample2\n")
            f.write("TP53,120.5,130.2\n")
            f.write("BRCA1,,50.1\n")  # 含 NA

        with open(cls.csv_corrupt, "wb") as f:
            f.write(b"\x00\x01\xFF\xFE\x00")  # 二进制数据，pd.read_csv 会报错

    @classmethod
    def tearDownClass(cls):
        import shutil
        shutil.rmtree(cls.tmpdir, ignore_errors=True)

    def _import_check(self):
        from GEO_expMatrix_check import check_single_matrix
        return check_single_matrix

    def test_clean_matrix(self):
        from GEO_expMatrix_check import check_single_matrix
        r = check_single_matrix(self.csv_clean, "GSE001")
        self.assertTrue(r["exists"])
        self.assertIsNone(r["error"])
        self.assertEqual(r["data_shape"], (2, 2))
        self.assertFalse(r["has_na"])
        self.assertFalse(r["has_negative"])
        self.assertAlmostEqual(r["min_value"], 45.0)
        self.assertAlmostEqual(r["max_value"], 130.2)

    def test_matrix_with_na(self):
        from GEO_expMatrix_check import check_single_matrix
        r = check_single_matrix(self.csv_with_na, "GSE002")
        self.assertTrue(r["exists"])
        self.assertTrue(r["has_na"])
        self.assertEqual(r["na_count"], 1)

    def test_corrupt_file(self):
        from GEO_expMatrix_check import check_single_matrix
        r = check_single_matrix(self.csv_corrupt, "GSE003")
        self.assertTrue(r["exists"])
        self.assertIsNotNone(r["error"])


class TestNCBI_TYPES(unittest.TestCase):
    """验证 NCBI_TYPES 常量。"""
    def test_types(self):
        self.assertIn("Count", NCBI_TYPES)
        self.assertIn("FPKM", NCBI_TYPES)
        self.assertIn("TPM", NCBI_TYPES)
        self.assertEqual(len(NCBI_TYPES), 3)


# ================================================================
# 从 GEO_smart_sniff.py 导入常量（不触发网络请求）
# ================================================================

# 手工加载常量，避免触发 module-level 的 Entrez 导入
def _load_sniff_constants():
    """安全加载 GEO_smart_sniff.py 中的纯常量，不执行网络代码。"""
    sniff_path = os.path.join(_SCRIPTS_DIR, "GEO_smart_sniff.py")
    with open(sniff_path, "r", encoding="utf-8") as f:
        source = f.read()

    # 提取常量定义（在 "Entrez.email" 赋值之前的部分，不会触发网络）
    consts = {}
    exec(compile(source, sniff_path, "exec"), consts)
    return consts


_SNIFF_CONSTS = _load_sniff_constants()
SPECIES_MAP          = _SNIFF_CONSTS.get("SPECIES_MAP", {})
MICROARRAY_PLATFORMS = _SNIFF_CONSTS.get("MICROARRAY_PLATFORMS", set())
RNA_SEQ_PLATFORMS    = _SNIFF_CONSTS.get("RNA_SEQ_PLATFORMS", set())
NCBI_GENOME_TAG      = _SNIFF_CONSTS.get("NCBI_GENOME_TAG", {})


class TestSpeciesMap(unittest.TestCase):
    """验证物种映射表完整性。"""

    def test_known_species(self):
        self.assertEqual(SPECIES_MAP["Homo sapiens"], ("human", "org.Hs.eg.db"))
        self.assertEqual(SPECIES_MAP["Mus musculus"], ("mouse", "org.Mm.eg.db"))
        self.assertEqual(SPECIES_MAP["Rattus norvegicus"], ("rat", "org.Rn.eg.db"))
        self.assertEqual(SPECIES_MAP["Danio rerio"], ("zebrafish", "org.Dr.eg.db"))

    def test_all_values_are_two_tuples(self):
        for latin, (short, orgdb) in SPECIES_MAP.items():
            self.assertIsInstance(latin, str)
            self.assertIsInstance(short, str)
            self.assertIsInstance(orgdb, str)
            self.assertTrue(orgdb.startswith("org."), f"bad OrgDb: {orgdb}")
            self.assertTrue(orgdb.endswith(".db"), f"bad OrgDb: {orgdb}")

    def test_no_overlap_between_platform_sets(self):
        overlap = MICROARRAY_PLATFORMS & RNA_SEQ_PLATFORMS
        self.assertEqual(overlap, set(), f"平台集合交叉: {overlap}")

    def test_genome_tags(self):
        self.assertEqual(NCBI_GENOME_TAG["human"], "GRCh38.p13")
        self.assertEqual(NCBI_GENOME_TAG["mouse"], "GRCm39")
        self.assertEqual(NCBI_GENOME_TAG["rat"], "mRatBN7.2")


# ================================================================
# GEO_run.py 测试
# ================================================================
class TestCheckDone(unittest.TestCase):
    """测试 GEO_run._check_done 函数。"""

    def setUp(self):
        self.tmpdir = tempfile.mkdtemp(prefix="test_run_")
        self.gse = "GSE99999"
        self.gse_dir = os.path.join(self.tmpdir, self.gse)
        os.makedirs(self.gse_dir)

    def tearDown(self):
        import shutil
        shutil.rmtree(self.tmpdir, ignore_errors=True)

    def _import(self):
        from GEO_run import _check_done
        return _check_done

    def test_has_expmatrix(self):
        """有 expMatrix.csv → True"""
        _check_done = self._import()
        open(os.path.join(self.gse_dir, "expMatrix.csv"), "w").close()
        self.assertTrue(_check_done(self.gse, self.tmpdir))

    def test_has_expmatrix_with_suffix(self):
        """有 expMatrix_GPL570.csv → True"""
        _check_done = self._import()
        open(os.path.join(self.gse_dir, "expMatrix_GPL570.csv"), "w").close()
        self.assertTrue(_check_done(self.gse, self.tmpdir))

    def test_no_expmatrix(self):
        """无 expMatrix → False"""
        _check_done = self._import()
        self.assertFalse(_check_done(self.gse, self.tmpdir))

    def test_other_files_only(self):
        """只有 clinical.csv 等非 expMatrix 文件 → False"""
        _check_done = self._import()
        open(os.path.join(self.gse_dir, "clinical.csv"), "w").close()
        open(os.path.join(self.gse_dir, "probe2gene.csv"), "w").close()
        self.assertFalse(_check_done(self.gse, self.tmpdir))


class TestSniffParse(unittest.TestCase):
    """测试 GEO_run._sniff 的 JSON 解析（mock subprocess）。"""

    def setUp(self):
        self.tmpdir = tempfile.mkdtemp(prefix="test_sniff_")

    def tearDown(self):
        import shutil
        shutil.rmtree(self.tmpdir, ignore_errors=True)

    def _import(self):
        from GEO_run import _sniff
        return _sniff

    @patch("GEO_run._run")
    def test_normal_json(self, mock_run):
        """正常 JSON 输出 → 正确解析"""
        mock_run.return_value = (True, json.dumps([
            {"gse_id": "GSE12345", "ok": True, "script": "GEO_download_probe.R",
             "args": ["--gse", "GSE12345"]}
        ]))
        _sniff = self._import()
        result = _sniff(["GSE12345"], "http://127.0.0.1:7897")
        self.assertIn("GSE12345", result)
        self.assertEqual(result["GSE12345"]["script"], "GEO_download_probe.R")

    @patch("GEO_run._run")
    def test_failed_gse_filtered(self, mock_run):
        """ok=False 的 GSE 被过滤"""
        mock_run.return_value = (True, json.dumps([
            {"gse_id": "GSE12345", "ok": False, "error": "未找到"},
            {"gse_id": "GSE67890", "ok": True, "script": "GEO_download_probe.R"}
        ]))
        _sniff = self._import()
        result = _sniff(["GSE12345", "GSE67890"], "http://127.0.0.1:7897")
        self.assertNotIn("GSE12345", result)
        self.assertIn("GSE67890", result)

    @patch("GEO_run._run")
    def test_invalid_json(self, mock_run):
        """无效 JSON 输出 → 返回空 dict"""
        mock_run.return_value = (True, "not a json at all")
        _sniff = self._import()
        result = _sniff(["GSE12345"], "http://127.0.0.1:7897")
        self.assertEqual(result, {})

    @patch("GEO_run._run")
    def test_subprocess_failure(self, mock_run):
        """subprocess 失败 → 返回空 dict"""
        mock_run.return_value = (False, "command not found")
        _sniff = self._import()
        result = _sniff(["GSE12345"], "http://127.0.0.1:7897")
        self.assertEqual(result, {})

    @patch("GEO_run._run")
    def test_empty_output(self, mock_run):
        """空输出 → 返回空 dict"""
        mock_run.return_value = (True, "")
        _sniff = self._import()
        result = _sniff(["GSE12345"], "http://127.0.0.1:7897")
        self.assertEqual(result, {})


class TestRunR(unittest.TestCase):
    """测试 GEO_run._run_r（mock subprocess + 产物检查）。"""

    def setUp(self):
        self.tmpdir = tempfile.mkdtemp(prefix="test_runr_")
        self.gse = "GSE99999"
        self.gse_dir = os.path.join(self.tmpdir, self.gse)
        os.makedirs(self.gse_dir)
        # 设置 ROOT（_run_r 内部用 str(ROOT)）
        import GEO_run
        GEO_run.ROOT = Path(self.tmpdir)

    def tearDown(self):
        import shutil
        shutil.rmtree(self.tmpdir, ignore_errors=True)

    def _import(self):
        from GEO_run import _run_r
        return _run_r

    def _make_info(self, script="GEO_download_probe.R", args=None):
        return {"script": script, "args": args or ["--gse", self.gse]}

    @patch("GEO_run._run")
    def test_success(self, mock_run):
        """R 退出码 0 + 产物存在 → True"""
        mock_run.return_value = (True, "完成")
        # 创建产物
        open(os.path.join(self.gse_dir, "expMatrix.csv"), "w").close()
        _run_r = self._import()
        result = _run_r(self.gse, self._make_info(), self.tmpdir, False, 600, "http://127.0.0.1:7897")
        self.assertTrue(result)

    @patch("GEO_run._run")
    def test_r_failure(self, mock_run):
        """R 退出码非 0 → False"""
        mock_run.return_value = (False, "Error: package not found")
        _run_r = self._import()
        result = _run_r(self.gse, self._make_info(), self.tmpdir, False, 600, "http://127.0.0.1:7897")
        self.assertFalse(result)

    @patch("GEO_run._run")
    def test_success_but_no_output(self, mock_run):
        """R 退出码 0 但无产物 → False"""
        mock_run.return_value = (True, "完成")
        # 不创建产物
        _run_r = self._import()
        result = _run_r(self.gse, self._make_info(), self.tmpdir, False, 600, "http://127.0.0.1:7897")
        self.assertFalse(result)

    @patch("GEO_run._run")
    def test_no_script(self, mock_run):
        """info 无 script → False（不调 subprocess）"""
        mock_run.return_value = (True, "")
        _run_r = self._import()
        result = _run_r(self.gse, {"script": "", "args": []}, self.tmpdir, False, 600, "")
        self.assertFalse(result)
        mock_run.assert_not_called()

    @patch("GEO_run._run")
    def test_diff_dedup(self, mock_run):
        """sniff args 已含 --diff TRUE 时不重复追加"""
        mock_run.return_value = (True, "ok")
        open(os.path.join(self.gse_dir, "expMatrix.csv"), "w").close()
        _run_r = self._import()
        info = self._make_info(args=["--gse", self.gse, "--diff", "TRUE"])
        _run_r(self.gse, info, self.tmpdir, True, 600, "")
        # 检查传给 _run 的命令中 --diff 只出现一次
        cmd = mock_run.call_args[0][0]
        diff_count = sum(1 for i, x in enumerate(cmd) if x == "--diff")
        self.assertEqual(diff_count, 1)


# ================================================================
# 运行
# ================================================================
if __name__ == "__main__":
    unittest.main(verbosity=2)
