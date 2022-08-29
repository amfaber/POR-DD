import asyncio

async def idk(cmd):
    proc = await asyncio.create_subprocess_shell(
        cmd,
        stdin = asyncio.subprocess.PIPE,
        stdout=asyncio.subprocess.PIPE,
        stderr=asyncio.subprocess.PIPE)

    async def wait_for_gnina():
        output = ""
        while output != b'--Chunk finished--\n':
            output = await proc.stdout.readline()
            # print(output)
    await wait_for_gnina()
    proc.stdin.write(b"Ready\n")
    await wait_for_gnina()
    proc.stdin.write(b"quit\n")
    await proc.wait()
    # return proc
cmd = "gnina -r /home/amfaber/Masters/POR-DD/data/raw_data/por_structures/3ES9_1.pdb -l /home/amfaber/Masters/POR-DD/data/raw_data/cyp_screen/first_50.sdf --minimize --continuous_operation -o gnina_out.sdf"
asyncio.run(idk(cmd))

