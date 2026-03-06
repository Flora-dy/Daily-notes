# Streamlit 部署说明

## 1. 本地启动

```bash
cd /Users/zzp/Code/Daily
streamlit run app.py
```

默认会打开 `http://localhost:8501`。

## 2. 推送到 GitHub

确保以下文件已提交到仓库：

- `app.py`
- `requirements.txt`
- `.streamlit/config.toml`
- `records.json`
- `ai-settings.json`

## 3. 部署到 Streamlit Community Cloud

1. 打开 [https://share.streamlit.io](https://share.streamlit.io)
2. 使用 GitHub 登录
3. 点击 **New app**
4. 选择你的仓库和分支
5. **Main file path** 填：`Daily/app.py`（如果仓库根目录就是 Daily，则填 `app.py`）
6. 点击 **Deploy**

## 4. 数据持久化说明（重要）

当前应用把数据写在：

- `records.json`
- `ai-settings.json`

在 Streamlit Community Cloud 上，容器重启后本地文件可能重置。  
如果要长期保存，建议下一步接入外部存储（如 Supabase / PostgreSQL / Google Sheets）。
