# 发布与推送步骤（GitHub → Zenodo DOI）

> ⚠️ **本文件为简版，已被 `outputs/TP投稿材料/09_GitHub与Zenodo发布操作手册.md` 取代**（含实测环境、5 项发布前清理、PAT 认证、Zenodo `.zenodo.json` 命名坑、完整回填清单）。**以 09 手册为准。**
>
> 一次性操作，约 10 分钟。完成后把 DOI 回填稿件 Declarations 与 cover letter。

## A. 推送 GitHub

```bash
cd repo_jakstat_mdd

# 1) 改 git 身份 + 【必须】重写已有的 7 个 commit 的作者
#    （当前全部为占位 jakstat-mdd <authors@example.com>；只改 config 不影响历史，
#     Zenodo 会从 GitHub contributor 统计取作者，占位会污染 DOI 页面）
git config user.name "你的真实姓名"
git config user.email "你的真实邮箱"
git filter-branch --env-filter '
  export GIT_AUTHOR_NAME="你的真实姓名"; export GIT_AUTHOR_EMAIL="你的真实邮箱"
  export GIT_COMMITTER_NAME="你的真实姓名"; export GIT_COMMITTER_EMAIL="你的真实邮箱"
' -- --all
git log --format='%an <%ae>' | sort | uniq -c   # 自检：不应再出现 authors@example.com

# 2) 在 github.com 新建空仓库（不要勾选 README/LICENSE），建议名：jakstat-mdd-gradient
#    认证：GitHub 已废除密码登录，需 Fine-grained PAT（Contents: read/write），见 09 手册 §3 步骤 4
git remote add origin https://github.com/xiwenHu/jakstat-mdd-gradient.git
git push -u origin master        # ← 本地分支是 master，不是 main
```

## B. 挂 Zenodo 拿 DOI

1. 打开 https://zenodo.org → 登录 → **GitHub** 标签 → 找到 `jakstat-mdd-gradient` → **On**（开启存档）
2. 回 GitHub：**Releases** → **Draft a new release** → tag `v1.0.0` → 发布
3. Zenodo 自动归档并生成 DOI（10.5281/zenodo.xxxxxxx）
4. 把 DOI 填回两处：
   - `zenodo.json` 的 `related_identifiers.identifier`（论文 DOI，接收后再填）
   - 本文件顶部 + README 的 Citation 段

## C. 稿件回填点（拿到 URL/DOI 后改这三处）

| 文件 | 位置 | 替换 |
|---|---|---|
| `docs/路径1_全文初稿.md` | Declarations | `[**补 GitHub URL**]` ×2 → 仓库 URL |
| `outputs/TP投稿材料/01_TP_cover_letter.md` | Declarations 段 | `[GitHub/Zenodo URL]` → Zenodo DOI 链接 |
| README | Citation 段 | BibTeX 的 `url`/`doi` 占位 |

## D. 备选：不想要 Zenodo

直接引用 GitHub 仓库也可（部分期刊接受），但 **Zenodo DOI 是期刊数据政策的硬通货**，TP 的 research data policy 倾向可引用的固定标识——建议走 B。
