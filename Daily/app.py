import io
import json
import re
import uuid
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd
import streamlit as st


ROOT_DIR = Path(__file__).resolve().parent
DATA_FILE = ROOT_DIR / "records.json"
SETTINGS_FILE = ROOT_DIR / "ai-settings.json"

PRIORITY_CN_TO_EN = {"普通": "normal", "高": "high", "紧急": "urgent"}
PRIORITY_EN_TO_CN = {"normal": "普通", "high": "高", "urgent": "紧急"}
STATUS_CN_TO_EN = {"待处理": "pending", "进行中": "in-progress", "已完成": "done"}
STATUS_EN_TO_CN = {"pending": "待处理", "in-progress": "进行中", "done": "已完成"}

DEFAULT_RULES = [
    {"keyword": "辅料", "category": "辅料策略"},
    {"keyword": "益生元", "category": "辅料策略"},
    {"keyword": "载体", "category": "辅料策略"},
    {"keyword": "功能辅料的数据库", "category": "数据库建设"},
    {"keyword": "数据库先积累", "category": "数据库建设"},
    {"keyword": "功能辅料", "category": "功能辅料研究"},
    {"keyword": "数据库", "category": "数据库建设"},
    {"keyword": "法规", "category": "法规评估"},
    {"keyword": "非转基因", "category": "质量合规"},
    {"keyword": "过敏原", "category": "质量合规"},
    {"keyword": "有机", "category": "质量合规"},
    {"keyword": "糖尿病", "category": "特殊人群"},
    {"keyword": "肿瘤患者", "category": "特殊人群"},
    {"keyword": "脑肠轴", "category": "产品开发"},
    {"keyword": "猴头菇", "category": "原料"},
    {"keyword": "临床", "category": "临床合作"},
    {"keyword": "配方", "category": "配方开发"},
    {"keyword": "周报", "category": "报告"},
    {"keyword": "月报", "category": "报告"},
    {"keyword": "汇报", "category": "报告"},
    {"keyword": "会议", "category": "会议"},
    {"keyword": "评审", "category": "会议"},
    {"keyword": "合同", "category": "法务"},
    {"keyword": "发票", "category": "财务"},
    {"keyword": "付款", "category": "财务"},
    {"keyword": "采购", "category": "采购"},
    {"keyword": "上线", "category": "项目"},
    {"keyword": "排期", "category": "项目"},
    {"keyword": "招聘", "category": "人事"},
]


def safe_text(value: str, max_len: int = 300) -> str:
    return str(value or "").strip()[:max_len]


def now_iso() -> str:
    return datetime.now(timezone.utc).isoformat().replace("+00:00", "Z")


def today_cn() -> str:
    now = datetime.now()
    return f"{now.month}月{now.day}日"


def normalize_category(value: str) -> str:
    raw = safe_text(value, 80)
    if not raw:
        return "通用"
    lower = raw.lower()
    if lower in {"general", "default", "other", "others", "misc", "normal"}:
        return "通用"
    if raw in {"通用", "默认", "其他", "未分类"}:
        return "通用"
    return raw


def normalize_source(value: str, fallback: str = "未标注") -> str:
    raw = safe_text(value, 80)
    if not raw:
        return fallback
    lower = raw.lower()
    if lower in {"ai", "assistant", "agent", "ai录入"}:
        return "AI录入"
    if lower in {"manual", "manual-input", "手动", "手工"}:
        return "手动录入"
    return raw


def normalize_priority(value: str) -> str:
    raw = safe_text(value, 20).lower()
    if raw in {"urgent", "u", "紧急", "特急", "高紧急", "asap", "red", "红色"}:
        return "urgent"
    if raw in {"high", "h", "高", "高优先", "高优先级", "important", "yellow"}:
        return "high"
    if raw in {"normal", "n", "普通", "一般", "low", "低"}:
        return "normal"
    return "normal"


def normalize_status(value: str) -> str:
    raw = safe_text(value, 30).lower()
    if raw in {"done", "completed", "finish", "finished", "已完成", "完成"}:
        return "done"
    if raw in {"in-progress", "in progress", "doing", "进行中", "处理中"}:
        return "in-progress"
    return "pending"


def read_records() -> list[dict]:
    if not DATA_FILE.exists():
        DATA_FILE.write_text("[]", encoding="utf-8")
        return []
    try:
        data = json.loads(DATA_FILE.read_text(encoding="utf-8"))
        return data if isinstance(data, list) else []
    except Exception:
        return []


def write_records(records: list[dict]) -> None:
    DATA_FILE.write_text(
        json.dumps(records, ensure_ascii=False, indent=2),
        encoding="utf-8",
    )


def read_rules() -> list[dict]:
    if not SETTINGS_FILE.exists():
        SETTINGS_FILE.write_text(
            json.dumps({"rules": DEFAULT_RULES}, ensure_ascii=False, indent=2),
            encoding="utf-8",
        )
        return DEFAULT_RULES
    try:
        data = json.loads(SETTINGS_FILE.read_text(encoding="utf-8"))
        rules = data.get("rules", [])
        return rules if isinstance(rules, list) and rules else DEFAULT_RULES
    except Exception:
        return DEFAULT_RULES


def write_rules(rules: list[dict]) -> None:
    SETTINGS_FILE.write_text(
        json.dumps({"rules": rules}, ensure_ascii=False, indent=2),
        encoding="utf-8",
    )


def infer_source_from_text(raw_text: str, fallback: str = "未标注") -> str:
    text = safe_text(raw_text, 1000)
    if not text:
        return fallback

    named_visit = re.search(
        r"([\u4e00-\u9fa5A-Za-z]{1,8}(?:总|经理|老师|博士)?)(?:拜访客户|走访客户|客户拜访)",
        text,
    )
    if named_visit:
        return f"{named_visit.group(1)}拜访客户"

    if re.search(r"(拜访客户|客户拜访|客户走访|走访客户)", text):
        return "客户拜访"
    if re.search(r"(客户反馈|客户提出|客户希望|客户要做|客户会自己做|客户沟通|客户需求|有客户)", text):
        return "客户沟通"
    if re.search(r"(领导|总经理|老板|安排|指示)", text):
        return "领导安排"
    if re.search(r"(会议纪要|会议|周会|例会|评审会)", text):
        return "会议纪要"
    if re.search(r"(邮件|邮件沟通)", text):
        return "邮件沟通"
    if re.search(r"(电话|电话沟通|通话)", text):
        return "电话沟通"
    return fallback


def display_source(source: str, demand: str, note: str) -> str:
    raw = safe_text(source, 80)
    if raw and raw not in {"未标注", "AI录入"}:
        return raw
    inferred = infer_source_from_text(f"{demand} {note}", fallback="未标注")
    if inferred in {"", "未标注", "AI录入"}:
        return "未分解"
    return inferred


def infer_priority_from_text(raw_text: str) -> str:
    text = safe_text(raw_text, 1000)
    if re.search(r"(紧急|马上|立刻|今天必须|当天|加急|asap)", text, re.I):
        return "urgent"
    if re.search(r"(尽快|优先|本周|抓紧|尽早)", text, re.I):
        return "high"
    return "normal"


def infer_status_from_text(raw_text: str) -> str:
    text = safe_text(raw_text, 1000)
    if re.search(r"(已完成|完成了|已提交|已发送|处理完|done|completed)", text, re.I):
        return "done"
    if re.search(r"(进行中|推进中|跟进中|在做|处理中|in progress|in-progress)", text, re.I):
        return "in-progress"
    if re.search(r"(继续|持续|先积累|继续挖|深化|推进)", text):
        return "in-progress"
    return "pending"


def infer_note_from_text(raw_text: str) -> str:
    notes = []
    if re.search(r"(别着急|不着急|暂缓)", raw_text):
        notes.append("该事项可暂缓推进")
    if re.search(r"(客户会自己做|客户自行处理|客户自行完成)", raw_text):
        notes.append("客户会自行处理部分内容")
    if re.search(r"(一起做临床|联合临床|可做临床|可以一起做临床)", raw_text):
        notes.append("可与客户联合开展临床")
    return "；".join(notes)


def strip_leading_date(text: str) -> str:
    value = safe_text(text, 1200)
    return re.sub(
        r"^\s*(\d{1,2}月\d{1,2}日|[12]\d{3}[年./-]\d{1,2}[月./-]\d{1,2}日?|[12]\d{3}/\d{1,2}/\d{1,2})\s*[:：]?\s*",
        "",
        value,
    )


def detect_builtin_category(text: str) -> str:
    mapping = [
        (["数据库", "积累数据", "数据沉淀", "数据建设"], "数据库建设"),
        (["辅料", "益生元", "载体"], "辅料策略"),
        (["功能辅料", "原料积累", "功能成分"], "功能辅料研究"),
        (["糖尿病", "肿瘤", "特殊人群"], "特殊人群"),
        (["法规", "过敏原", "非转基因", "有机", "地域特色"], "质量合规"),
        (["脑肠轴", "猴头菇", "记忆力", "配方"], "产品开发"),
        (["临床", "受试", "试验"], "临床合作"),
        (["客户拜访", "客户反馈", "客户需求"], "客户需求"),
        (["周报", "月报", "汇报", "报告", "PPT"], "报告"),
        (["会议", "评审", "例会", "会审"], "会议"),
        (["合同", "法务", "条款"], "法务"),
        (["发票", "报销", "付款", "对账"], "财务"),
        (["采购", "下单", "供应商"], "采购"),
        (["上线", "需求", "排期", "里程碑", "版本"], "项目"),
        (["招聘", "面试", "入职", "离职"], "人事"),
    ]
    for keywords, category in mapping:
        if any(keyword in text for keyword in keywords):
            return category
    return "通用"


def detect_category_by_rules(text: str, rules: list[dict]) -> str:
    matched = ""
    max_len = 0
    for rule in rules:
        keyword = safe_text(rule.get("keyword", ""), 40)
        category = normalize_category(rule.get("category", ""))
        if not keyword or not category:
            continue
        if keyword in text and len(keyword) > max_len:
            matched = category
            max_len = len(keyword)
    return matched


def split_ai_input_to_items(text: str) -> list[str]:
    normalized = (
        str(text or "")
        .replace("\r\n", "\n")
        .replace("\u3000", " ")
    )
    normalized = "\n".join(
        re.sub(r"[ \t]+", " ", line).strip() for line in normalized.split("\n")
    ).strip()
    if not normalized:
        return []

    lines_raw = [line.strip() for line in normalized.split("\n") if line.strip()]
    global_source = ""
    if lines_raw and re.match(r"^(?:来源|source)\s*[:：]", lines_raw[0], re.I):
        m = re.match(r"^(?:来源|source)\s*[:：]?\s*(.+)$", lines_raw[0], re.I)
        if m:
            global_source = safe_text(m.group(1), 40)
    lines = [
        line
        for idx, line in enumerate(lines_raw)
        if not (idx == 0 and re.match(r"^(?:来源|source)\s*[:：]", line, re.I))
    ]

    items = []
    current = ""
    for line in lines:
        if re.match(r"^(\d+[\.\)、．]|[-*•])\s*", line):
            if current:
                items.append(current.strip())
            current = re.sub(r"^(\d+[\.\)、．]|[-*•])\s*", "", line).strip()
        else:
            current = f"{current} {line}".strip() if current else line
    if current:
        items.append(current.strip())

    def apply_global_source(input_items: list[str]) -> list[str]:
        if not global_source:
            return input_items
        result = []
        for item in input_items:
            if re.search(r"(?:来源|source)\s*[:：]", item, re.I):
                result.append(item)
            else:
                result.append(f"{item} 来源:{global_source}")
        return result

    if len(items) >= 2:
        return apply_global_source([item for item in items if item])

    markers = list(re.finditer(r"(^|[^\d])(\d+[\.\)、．]\s*)", normalized))
    inline = []
    if len(markers) >= 2:
        for idx, marker in enumerate(markers):
            marker_start = marker.start() + len(marker.group(1))
            content_start = marker_start + len(marker.group(2))
            content_end = (
                markers[idx + 1].start() + len(markers[idx + 1].group(1))
                if idx + 1 < len(markers)
                else len(normalized)
            )
            piece = normalized[content_start:content_end].strip()
            if piece:
                inline.append(piece)
    if len(inline) >= 2:
        return apply_global_source(inline)

    return apply_global_source([normalized])


def parse_ai_record(text: str, rules: list[dict]) -> dict:
    raw = safe_text(text, 1200)
    if not raw:
        return {}

    explicit_source = ""
    source_match = re.search(r"(?:来源|source)\s*[:：]?\s*([^\n]{1,120})", raw, re.I)
    if source_match:
        explicit_source = safe_text(
            re.sub(
                r"^\s*(\d{1,2}月\d{1,2}日|[12]\d{3}[年./-]\d{1,2}[月./-]\d{1,2}日?)\s*[:：]?\s*",
                "",
                source_match.group(1),
            ),
            80,
        )
    source = (
        normalize_source(explicit_source, "AI录入")
        if explicit_source
        else infer_source_from_text(raw, "AI录入")
    )

    priority = infer_priority_from_text(raw)
    priority_match = re.search(r"(?:优先级|priority)\s*[:：]\s*([\u4e00-\u9fa5A-Za-z0-9_-]{1,20})", raw, re.I)
    if priority_match:
        priority = normalize_priority(priority_match.group(1))

    status = infer_status_from_text(raw)
    status_match = re.search(r"(?:状态|status)\s*[:：]\s*([\u4e00-\u9fa5A-Za-z0-9_-]{1,20})", raw, re.I)
    if status_match:
        status = normalize_status(status_match.group(1))

    category = ""
    category_match = re.search(r"(?:分类|category)\s*[:：]\s*([\u4e00-\u9fa5A-Za-z0-9_-]{1,30})", raw, re.I)
    if category_match:
        category = normalize_category(category_match.group(1))
    if not category or category == "通用":
        category = detect_category_by_rules(raw, rules) or detect_builtin_category(raw)

    note = ""
    note_match = re.search(r"(?:备注|note)\s*[:：]\s*(.+)$", raw, re.I)
    if note_match:
        note = safe_text(note_match.group(1), 500)
    else:
        note = safe_text(infer_note_from_text(raw), 500)

    demand = strip_leading_date(raw)
    demand = re.sub(r"^\s*\d+[\.\)、．]\s*", "", demand)
    demand = re.sub(r"(?:来源|source)\s*[:：]\s*[^\n，,。;；]{1,40}", "", demand, flags=re.I)
    demand = re.sub(r"(?:分类|category)\s*[:：]\s*[\u4e00-\u9fa5A-Za-z0-9_-]{1,30}", "", demand, flags=re.I)
    demand = re.sub(r"(?:优先级|priority)\s*[:：]\s*[\u4e00-\u9fa5A-Za-z0-9_-]{1,30}", "", demand, flags=re.I)
    demand = re.sub(r"(?:状态|status)\s*[:：]\s*[\u4e00-\u9fa5A-Za-z0-9_-]{1,30}", "", demand, flags=re.I)
    demand = re.sub(r"(?:备注|note)\s*[:：].+$", "", demand, flags=re.I)
    demand = re.sub(
        r"^[\u4e00-\u9fa5A-Za-z]{2,12}(总|经理|老师|博士)?(?:拜访客户|反馈|提出|沟通|联系)\s*[，,:：]*",
        "",
        demand,
        flags=re.I,
    ).strip()
    demand_match = re.search(r"(?:客户要做|客户提出|客户希望|需要|需做|请)(.+)", demand)
    if demand_match:
        demand = f"跟进{demand_match.group(1).strip()}"
    demand = re.sub(r"^[，,。;；\s]+", "", demand)
    demand = re.sub(r"[，,。;；\s]+$", "", demand)
    demand = re.sub(r"\s{2,}", " ", demand)
    if not demand:
        demand = f"待补充任务（{today_cn()}）"

    return {
        "demand": safe_text(demand, 300),
        "category": normalize_category(category),
        "priority": normalize_priority(priority),
        "status": normalize_status(status),
        "source": normalize_source(source, "AI录入"),
        "note": safe_text(note, 500),
    }


def build_record(input_data: dict) -> dict:
    now = now_iso()
    return {
        "id": f"{int(datetime.now().timestamp() * 1000)}-{uuid.uuid4().hex[:5]}",
        "demand": safe_text(input_data.get("demand", ""), 300),
        "category": normalize_category(input_data.get("category", "")),
        "priority": normalize_priority(input_data.get("priority", "normal")),
        "status": normalize_status(input_data.get("status", "pending")),
        "source": normalize_source(input_data.get("source", ""), "手动录入"),
        "note": safe_text(input_data.get("note", ""), 500),
        "createdAt": now,
        "updatedAt": now,
    }


def rules_to_text(rules: list[dict]) -> str:
    return "\n".join(
        f"{safe_text(item.get('keyword', ''), 40)}={safe_text(item.get('category', ''), 40)}"
        for item in rules
        if safe_text(item.get("keyword", ""), 40) and safe_text(item.get("category", ""), 40)
    )


def parse_rules_text(text: str) -> list[dict]:
    rules = []
    for line in str(text or "").splitlines():
        line = line.strip()
        if not line:
            continue
        match = re.match(r"^(.+?)[=:：](.+)$", line)
        if not match:
            continue
        keyword = safe_text(match.group(1), 40)
        category = normalize_category(match.group(2))
        if keyword and category:
            rules.append({"keyword": keyword, "category": category})
    return rules


def to_excel_bytes(records: list[dict]) -> bytes:
    rows = []
    for idx, rec in enumerate(records, start=1):
        rows.append(
            {
                "序号": idx,
                "事项": rec.get("demand", ""),
                "分类": rec.get("category", ""),
                "优先级": PRIORITY_EN_TO_CN.get(rec.get("priority", "normal"), "普通"),
                "来源": display_source(rec.get("source", ""), rec.get("demand", ""), rec.get("note", "")),
                "状态": STATUS_EN_TO_CN.get(rec.get("status", "pending"), "待处理"),
                "备注": rec.get("note", ""),
                "创建时间": rec.get("createdAt", ""),
                "更新时间": rec.get("updatedAt", ""),
            }
        )
    df = pd.DataFrame(rows)
    output = io.BytesIO()
    with pd.ExcelWriter(output, engine="openpyxl") as writer:
        df.to_excel(writer, sheet_name="流水记录", index=False)
    output.seek(0)
    return output.getvalue()


def filter_records(records: list[dict], status_filter: str, priority_filter: str, keyword: str) -> list[dict]:
    filtered = records
    if status_filter != "全部":
        filtered = [r for r in filtered if STATUS_EN_TO_CN.get(r.get("status", "pending"), "待处理") == status_filter]
    if priority_filter != "全部":
        filtered = [r for r in filtered if PRIORITY_EN_TO_CN.get(r.get("priority", "normal"), "普通") == priority_filter]
    key = safe_text(keyword, 80)
    if key:
        filtered = [
            r
            for r in filtered
            if key in f"{r.get('demand', '')} {r.get('category', '')} {r.get('source', '')} {r.get('note', '')}"
        ]
    return filtered


st.set_page_config(page_title="日常流水记录助手", layout="wide")
st.markdown(
    """
<style>
.stApp { background: #f4f7f2; }
div[data-testid="stMetric"] { background: #ffffff; border: 1px solid #dfe7d9; padding: 8px 10px; }
</style>
""",
    unsafe_allow_html=True,
)

st.title("日常流水记录助手（Streamlit）")
st.caption("支持手动登记、AI 批量拆分、实时编辑、导出 Excel")

records = read_records()
rules = read_rules()

with st.sidebar:
    st.subheader("筛选")
    status_filter = st.selectbox("状态", ["全部", "待处理", "进行中", "已完成"], index=0)
    priority_filter = st.selectbox("优先级", ["全部", "普通", "高", "紧急"], index=0)
    keyword = st.text_input("关键词搜索")

    st.divider()
    st.subheader("AI 自动分类规则")
    rules_text = st.text_area(
        "格式：关键词=分类（每行一条）",
        value=rules_to_text(rules),
        height=260,
    )
    if st.button("保存分类规则", use_container_width=True):
        new_rules = parse_rules_text(rules_text)
        if not new_rules:
            st.error("至少填写一条有效规则，例如：周报=报告")
        else:
            write_rules(new_rules)
            st.success(f"已保存 {len(new_rules)} 条规则")
            st.rerun()

total = len(records)
urgent = sum(1 for r in records if r.get("priority") == "urgent")
pending = sum(1 for r in records if r.get("status") == "pending")
done = sum(1 for r in records if r.get("status") == "done")

c1, c2, c3, c4 = st.columns(4)
c1.metric("总数", total)
c2.metric("紧急", urgent)
c3.metric("待处理", pending)
c4.metric("已完成", done)

tab_add, tab_list, tab_edit = st.tabs(["新增记录", "流水列表", "编辑记录"])

with tab_add:
    with st.form("manual_add", clear_on_submit=True):
        st.subheader("手动新增")
        demand = st.text_input("事项（必填）")
        col_a, col_b = st.columns(2)
        category = col_a.text_input("分类（可选）")
        source = col_b.text_input("来源（可选）")
        note = st.text_input("备注（可选）")
        col_c, col_d = st.columns(2)
        priority_cn = col_c.selectbox("优先级", ["普通", "高", "紧急"], index=0)
        status_cn = col_d.selectbox("状态", ["待处理", "进行中", "已完成"], index=0)
        submit_manual = st.form_submit_button("新增记录", use_container_width=True)
        if submit_manual:
            if not safe_text(demand, 300):
                st.error("事项不能为空")
            else:
                new_record = build_record(
                    {
                        "demand": demand,
                        "category": category,
                        "source": source,
                        "note": note,
                        "priority": PRIORITY_CN_TO_EN[priority_cn],
                        "status": STATUS_CN_TO_EN[status_cn],
                    }
                )
                latest = read_records()
                latest.insert(0, new_record)
                write_records(latest)
                st.success("记录已新增")
                st.rerun()

    st.divider()
    with st.form("ai_add", clear_on_submit=True):
        st.subheader("AI 新增")
        ai_text = st.text_area("可粘贴一条或多条（支持 1. 2. 3. 自动拆分）", height=180)
        submit_ai = st.form_submit_button("AI 拆分并入库", use_container_width=True)
        if submit_ai:
            text = safe_text(ai_text, 4000)
            if not text:
                st.error("请先输入要解析的文本")
            else:
                chunks = split_ai_input_to_items(text)[:20]
                created = []
                for chunk in chunks:
                    parsed = parse_ai_record(chunk, rules)
                    if parsed.get("demand"):
                        created.append(build_record(parsed))
                if not created:
                    st.error("未解析出有效事项，请调整文本后重试")
                else:
                    latest = read_records()
                    write_records(created + latest)
                    st.success(f"已拆分 {len(chunks)} 段，新增 {len(created)} 条记录")
                    st.rerun()

with tab_list:
    filtered_records = filter_records(records, status_filter, priority_filter, keyword)
    st.subheader(f"流水列表（{len(filtered_records)} 条）")
    display_rows = []
    for rec in filtered_records:
        display_rows.append(
            {
                "时间": safe_text(rec.get("createdAt", ""), 30).replace("T", " ")[:16],
                "事项": rec.get("demand", ""),
                "来源": display_source(rec.get("source", ""), rec.get("demand", ""), rec.get("note", "")),
                "分类": rec.get("category", ""),
                "优先级": PRIORITY_EN_TO_CN.get(rec.get("priority", "normal"), "普通"),
                "状态": STATUS_EN_TO_CN.get(rec.get("status", "pending"), "待处理"),
                "备注": rec.get("note", ""),
                "ID": rec.get("id", ""),
            }
        )
    st.dataframe(pd.DataFrame(display_rows), use_container_width=True, hide_index=True)

    excel_data = to_excel_bytes(filtered_records)
    st.download_button(
        "导出当前列表为 Excel",
        data=excel_data,
        file_name=f"daily-records-{datetime.now().strftime('%Y-%m-%d')}.xlsx",
        mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet",
        use_container_width=True,
    )

with tab_edit:
    st.subheader("编辑 / 删除")
    if not records:
        st.info("暂无记录可编辑")
    else:
        options = {
            f"{safe_text(rec.get('createdAt', ''), 30).replace('T', ' ')[:16]} | {safe_text(rec.get('demand', ''), 40)}": rec["id"]
            for rec in records
        }
        selected_label = st.selectbox("选择记录", list(options.keys()))
        selected_id = options[selected_label]
        target_idx = next((i for i, rec in enumerate(records) if rec.get("id") == selected_id), -1)
        target = records[target_idx] if target_idx >= 0 else None

        if target:
            default_source = safe_text(target.get("source", ""), 80)
            if default_source in {"AI录入", "未标注"}:
                default_source = ""

            with st.form("edit_form"):
                edit_demand = st.text_input("事项", value=target.get("demand", ""))
                col_e, col_f = st.columns(2)
                edit_category = col_e.text_input("分类", value=target.get("category", ""))
                edit_source = col_f.text_input("来源", value=default_source)
                edit_note = st.text_input("备注", value=target.get("note", ""))
                col_g, col_h = st.columns(2)
                edit_priority_cn = col_g.selectbox(
                    "优先级",
                    ["普通", "高", "紧急"],
                    index=["normal", "high", "urgent"].index(target.get("priority", "normal")),
                )
                edit_status_cn = col_h.selectbox(
                    "状态",
                    ["待处理", "进行中", "已完成"],
                    index=["pending", "in-progress", "done"].index(target.get("status", "pending")),
                )
                save_edit = st.form_submit_button("保存修改", use_container_width=True)
                if save_edit:
                    demand_val = safe_text(edit_demand, 300)
                    if not demand_val:
                        st.error("事项不能为空")
                    else:
                        records[target_idx]["demand"] = demand_val
                        records[target_idx]["category"] = normalize_category(edit_category)
                        records[target_idx]["priority"] = PRIORITY_CN_TO_EN[edit_priority_cn]
                        records[target_idx]["status"] = STATUS_CN_TO_EN[edit_status_cn]
                        src_raw = safe_text(edit_source, 80)
                        if not src_raw:
                            src_raw = infer_source_from_text(f"{edit_demand} {edit_note}", fallback="未标注")
                        records[target_idx]["source"] = normalize_source(src_raw, "未标注")
                        records[target_idx]["note"] = safe_text(edit_note, 500)
                        records[target_idx]["updatedAt"] = now_iso()
                        write_records(records)
                        st.success("保存成功")
                        st.rerun()

            if st.button("删除这条记录", type="secondary", use_container_width=True):
                latest = [rec for rec in records if rec.get("id") != selected_id]
                write_records(latest)
                st.success("已删除")
                st.rerun()
