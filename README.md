# Paired-end Fragment QC Review Service

建库质控复核服务：对一个**成功配对任务**（任务号随请求提交）的双端读对做
异常片段复核。输入两份配对 FASTQ（R1/R2）与覆盖**全部读对**的比对候选
JSON；服务先校准正常片段跨度，再把每对判为正常、跨参考、朝向异常、跨度
异常或待定，输出异常读号、类别及各类计数，供结构变异复核。FastAPI 实现，
单进程内存上限 256 MiB，无状态、不落盘、不改原件。

## 输入约束

- 至多 **2000 对**；每端序列不超过 **500 碱基**。
- 每个读段的每端至多 **3 个**候选；候选字段：
  `read_id`、`side`（1=R1，2=R2）、`reference`（参考名称）、
  `start`（**1 基起点**，正整数）、`strand`（`+` 正向 / `-` 反向）、
  `score`（整数得分）。
- 比对按**无插缺的全长比对**处理，即末端位置 = 起点 + 该端读长 − 1。
- 候选 JSON 形状：

```json
{
  "candidates": [
    {"read_id": "r0001", "side": 1, "reference": "chr1",
     "start": 1000, "strand": "+", "score": 50},
    {"read_id": "r0001", "side": 2, "reference": "chr1",
     "start": 1200, "strand": "-", "score": 48}
  ]
}
```

## 校准规则

仅从同时满足以下条件的「干净对子」中取样：

1. 两端各只有 1 个候选；
2. 两个候选位于**同一参考**；
3. R1 为正向且在左、R2 为反向且在右（朝内：R1 `+`、R2 `-`、
   `start_R1 < start_R2`）。

- 干净对子 **不足 20 对 → 整次请求报校准失败**（HTTP 400，不给分类结果）。
- 跨度定义（朝内全长比对）：
  `span = (start_R2 + len_R2 − 1) − start_R1 + 1`。
- 计算跨度的**中位数**与 **MAD**（中位数绝对偏差）；
  容差 `tolerance = max(50, 4 × MAD)`。
- 朝内对子满足 `|span − 中位数| ≤ tolerance` 判正常，否则判**跨度异常**
  （边界等值算正常）。

## 分类与最高分仲裁

按以下顺序给每个候选组合下结论：

| 情形 | 类别 |
| --- | --- |
| 两端参考不同 | `cross_reference`（跨参考） |
| 同参考但不朝内（链向或左右顺序不符） | `orientation`（朝向异常） |
| 朝内且跨度在容差内 | `normal` |
| 朝内但跨度超出容差 | `span`（跨度异常） |

- 每对枚举 R1×R2 全部候选组合，取**两端得分之和最高**的组合。
- 最高分并列时：仅当所有并列组合的结论（normal / cross_reference /
  orientation / span）**完全一致**才采信该结论；否则记为
  `undetermined`（待定）。

## 整次失败，不给部分结果

以下任一情况，整个请求返回错误（400/422），响应里不含任何分类结果，
且上传原件只读、绝不改写：

- 候选**漏号**（某读某端无候选）、引用了配对文件中不存在的读号、
  同一身份候选**重复**、**起点非正**（< 1）、每端候选超过 3 个；
- 候选 JSON 非法、字段缺失或类型错误（`side` 非 1/2、链向非 `+/-`、
  得分非整数等）；
- 配对 FASTQ **缺件/空文件**、非四行式、头行不以 `@` 开头、第三行不以
  `+` 开头、序列与质量不等长、**后段截断损坏**；
- 两侧读号集合/顺序不一致、侧内读号重复；
- 读对数超过 2000、任一读端超过 500 碱基；
- 校准对子不足 20。

## API

`POST /qc/analyze`，`multipart/form-data`：

| 字段 | 说明 |
| --- | --- |
| `task_id` | 成功配对任务的任务号（表单字段） |
| `r1` | 配对后的 R1 FASTQ（未压缩） |
| `r2` | 配对后的 R2 FASTQ（未压缩） |
| `candidates` | 候选 JSON 文件 |

响应 `200`：

```json
{
  "pair_count": 28,
  "calibration": {
    "pair_count": 24,
    "median_span": 300.0,
    "mad_span": 0.0,
    "span_tolerance": 50.0
  },
  "counts": {
    "normal": 24,
    "cross_reference": 1,
    "orientation": 1,
    "span": 1,
    "undetermined": 1
  },
  "anomalies": [
    {"position": 25, "read_id": "xref", "category": "cross_reference"},
    {"position": 26, "read_id": "orient", "category": "orientation"},
    {"position": 27, "read_id": "spanout", "category": "span"},
    {"position": 28, "read_id": "tied", "category": "undetermined"}
  ],
  "task_id": "demo-paired-task-001"
}
```

错误响应（无任何部分结果）：`{"error": "...", "detail": "..."}`。

另有 `GET /health` 健康检查；启动后打开交互式文档
<http://localhost:8000/docs> 可直接在页面上提交请求。

## 运行（Docker Compose）

```bash
docker compose up --build
# 打开 http://localhost:8000/docs
```

`docker-compose.yml` 中 `mem_limit: 256m` 为进程内存上限；服务无状态，
不挂载数据卷，请求处理完即释放全部资源。

## 本地直接运行（开发/自测）

```bash
python3 -m venv --without-pip .venv
.venv/bin/python tools/get-pip.py
.venv/bin/pip install -r requirements.txt
.venv/bin/uvicorn app.main:app --port 8000
```

## 演示脚本

```bash
./scripts/demo.sh                              # 需要服务运行在 localhost:8000
BASE_URL=http://localhost:8000 ./scripts/demo.sh
```

脚本自动生成配对 FASTQ 与候选文件，演示一次含四类异常的完整分析，以及
配对文件损坏、起点非正、校准不足三种整次失败。

## 自测

```bash
bash scripts/run_tests.sh
# 或
.venv/bin/python -I tests/test_qc.py
```

共 **60 项检查**，覆盖：中位数/MAD/容差校准与 20 对下限、跨度边界判定、
跨参考/朝向/跨度分类、最高分仲裁、并列一致采信与冲突待定、全部输入校验
错误、2000 对/500 碱基/3 候选上限、FASTQ 解析与配对校验，以及 HTTP 全链路
（200/400/422、错误信封无部分结果、`/health`、`/docs`、OpenAPI 注册）。

## 目录结构

```
app/
  config.py    # 2000 对 / 500 碱基 / 3 候选 / 20 校准对等上限（环境变量可覆盖）
  fastq.py     # FASTQ 解析、四行式校验、双端配对校验（只读原件）
  qc.py        # 候选解析与校验、校准、最高分仲裁、分类（纯逻辑）
  main.py      # FastAPI 路由、统一错误信封、有界上传
tests/test_qc.py
scripts/demo.sh, scripts/run_tests.sh
Dockerfile, docker-compose.yml, requirements.txt
```

## 边界说明 / 遗留项

- 仅接受未压缩 FASTQ；`.gz` 等压缩格式不在范围内。
- 比对固定按无插缺全长比对处理，不接受 CIGAR；片段跨度即两端比对覆盖的
  参考碱基数。
- 服务无持久化设计，请求结果不保存；如需留存复核结论，由调用方自行归档。
- 开发沙箱内未安装 Docker，`docker compose up` 未实跑；compose/Dockerfile
  仅做了 YAML 与结构校验，应用本身已用等价的本地 uvicorn 命令完整验证
  （含 `/docs`、200/400/422 全部分支）。在装有 Docker 的环境执行
  `docker compose up --build` 即可。
