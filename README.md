# PCR Primer Dimer Screening API

筛查 PCR 引物二聚体的 REST API。调用方提交一批带名称的 IUPAC DNA 引物及两个阈值
（连续互补长度、3′ 端互补长度），服务枚举每条引物的**自二聚体**和全部两两
**异二聚体**的所有错位对齐，按 IUPAC 简并碱基集合判断互补，对每对引物返回**最危险
对齐**及其风险原因。所有批次结果持久化在 SQLite 中，可随时查询。

## 快速启动（单条命令）

```bash
docker compose up --build
```

启动后即可访问：

| 地址 | 说明 |
| --- | --- |
| http://localhost:8000 | API 服务（端口 **8000**） |
| http://localhost:8000/docs | Swagger 交互式文档 |
| http://localhost:8000/health | 健康检查 |

依赖安装与数据库初始化均已纳入启动流程：容器构建时安装 Python 依赖，应用启动时
自动创建 SQLite 数据库文件与表结构，无需任何手工步骤。

## 配置

通过环境变量配置（`docker-compose.yml` 中已给出默认值）：

| 变量 | 默认值 | 说明 |
| --- | --- | --- |
| `DATABASE_PATH` | `/data/primer_dimer.db`（容器）/ `./primer_dimer.db`（本地） | SQLite 数据库文件路径，目录不存在时自动创建 |

端口映射在 `docker-compose.yml` 的 `ports` 中修改（默认 `8000:8000`）。数据库文件
通过命名卷 `primer-data` 持久化，容器重建后批次数据不丢失。

## API 使用

### 提交筛查批次

```bash
curl -X POST http://localhost:8000/api/v1/batches \
  -H 'Content-Type: application/json' \
  -d '{
    "primers": [
      {"name": "fwd", "sequence": "ATGCATGC"},
      {"name": "rev", "sequence": "GCATGCAT"}
    ],
    "max_consecutive_threshold": 4,
    "three_prime_threshold": 3
  }'
```

响应（`201 Created`）：

```json
{
  "batch_id": 1,
  "created_at": "2026-09-18T20:45:21.933392+00:00",
  "max_consecutive_threshold": 4,
  "three_prime_threshold": 3,
  "primer_count": 2,
  "pair_count": 3,
  "flagged_count": 3,
  "primers": [{"name": "fwd", "sequence": "ATGCATGC"}, ...],
  "results": [
    {
      "primer_a": "fwd",
      "primer_b": "rev",
      "kind": "hetero",
      "offset": 0,
      "max_consecutive": 8,
      "three_prime_consecutive": 8,
      "flagged": true,
      "reasons": ["MAX_CONSECUTIVE_COMPLEMENTARITY", "THREE_PRIME_END_COMPLEMENTARITY"]
    }
  ]
}
```

### 查询批次结果

```bash
curl http://localhost:8000/api/v1/batches/1      # 单个批次的完整报告
curl http://localhost:8000/api/v1/batches        # 所有批次摘要（新的在前）
```

### 端点一览

| 方法 | 路径 | 说明 |
| --- | --- | --- |
| `POST` | `/api/v1/batches` | 提交一批引物，同步筛查、入库并返回完整报告 |
| `GET` | `/api/v1/batches/{batch_id}` | 查询某批次的完整筛查报告（不存在返回 404） |
| `GET` | `/api/v1/batches` | 列出全部批次摘要 |
| `GET` | `/health` | 健康检查 |

### 请求参数

| 字段 | 类型 | 默认 | 约束 |
| --- | --- | --- | --- |
| `primers[].name` | string | 必填 | 1–128 字符，批次内唯一 |
| `primers[].sequence` | string | 必填 | 1–256 字符，仅允许 IUPAC 符号（大小写均可，自动归一为大写） |
| `primers` | array | 必填 | 1–500 条 |
| `max_consecutive_threshold` | int | 4 | ≥1，连续互补长度达到该值即标记风险 |
| `three_prime_threshold` | int | 3 | ≥1，3′ 端连续互补长度达到该值即标记风险 |

非法碱基符号、重名、空批次等返回 `422` 及具体错误信息。

## 筛查算法

**IUPAC 互补判定。** 每个简并碱基视为其代表的具体碱基集合（如 `R={A,G}`、
`N={A,C,G,T}`、`U` 按 `T` 处理）。两个碱基互补，当且仅当两个集合中各存在至少
一个可形成 Watson-Crick 配对（A-T 或 C-G）的具体碱基。

**对齐枚举。** 二聚体为反向平行结合：将引物 B 反向（`rev_b`，其索引 0 即 B 的
3′ 端）后与引物 A 做错位对齐。偏移量 `offset` 表示与 `seq_a[0]` 配对的
`rev_b` 索引，枚举范围为 `-(len(B)-1)` 到 `len(A)-1)`，覆盖全部错位。自二聚体
只枚举 `offset >= 0`：两条相同分子互换后偏移 `o` 与 `-o` 是同一物理二聚体，避免
重复计数。

**每对引物的指标。** 对每个对齐统计：

- `max_consecutive`：该对齐中最长的连续互补碱基段长度；
- `three_prime_consecutive`：涉及任一引物 3′ 端（A 的末端或 B 的 `rev_b[0]`，
  且该端落在重叠区内）的最长连续互补段长度。

**最危险对齐的选择（唯一结果）。** 先按 `max_consecutive` 最大者胜出；并列时按
`three_prime_consecutive` 最大者胜出；仍并列时按**固定偏移规则——偏移量最小者
胜出**，保证每对引物只报告唯一确定的对齐。

**风险原因（`reasons`）。**

| 代码 | 含义 |
| --- | --- |
| `MAX_CONSECUTIVE_COMPLEMENTARITY` | 最危险对齐的连续互补长度 ≥ `max_consecutive_threshold` |
| `THREE_PRIME_END_COMPLEMENTARITY` | 涉及 3′ 端的连续互补长度 ≥ `three_prime_threshold` |

`reasons` 非空时 `flagged=true`。结果按提交顺序确定性地排列：每条引物先报告自
二聚体，再报告与后续引物的异二聚体。

## 本地开发

```bash
pip install -r requirements-dev.txt
python -m pytest tests/                 # 运行测试
uvicorn app.main:app --reload --port 8000   # 本地启动（数据库默认为 ./primer_dimer.db）
```

## 项目结构

```
app/
  main.py      FastAPI 应用与路由（启动时初始化数据库）
  dimer.py     核心算法：IUPAC 互补、对齐枚举、最危险对齐选择
  db.py        SQLite 连接与建表
  schemas.py   Pydantic 请求/响应模型与校验
tests/
  test_dimer.py  算法单元测试
  test_api.py    API 集成测试（TestClient + 临时数据库）
Dockerfile
docker-compose.yml
requirements.txt / requirements-dev.txt
```
