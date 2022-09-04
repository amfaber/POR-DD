import asyncio


class Waiter:

    @classmethod
    async def create(cls, cmd):
        self = cls()
        self.proc = await asyncio.create_subprocess_shell(
            cmd,
            stdin = asyncio.subprocess.PIPE,
            stdout=asyncio.subprocess.PIPE,
            stderr=asyncio.subprocess.PIPE)
        return self
    
    async def wait_for_gnina(self):
        output = ""
        i = 0
        while output != b'--Chunk finished--\n':
            output = await self.proc.stdout.readline()
            print(output)
            i += 1
            if i > 1000:
                break
    
    def signal_ready(self):
        self.stdin.write(b"Ready\n")
    
    async def shut_down_gnina(self):
        self.proc.stdin.write(b'quit\n')
        await self.proc.wait()
    

# async def idk(cmd):
#     proc = await asyncio.create_subprocess_shell(
#         cmd,
#         stdin = asyncio.subprocess.PIPE,
#         stdout=asyncio.subprocess.PIPE,
#         stderr=asyncio.subprocess.PIPE)

#     async def wait_for_gnina():
#         output = ""
#         while output != b'--Chunk finished--\n':
#             output = await proc.stdout.readline()
#             # print(output)
#     await wait_for_gnina()
#     proc.stdin.write(b"Ready\n")
#     await wait_for_gnina()
#     proc.stdin.write(b"quit\n")
#     await proc.wait()
#     return proc

async def idk(cmd):
    waiter = await Waiter.create(cmd)
    await waiter.wait_for_gnina()
    await waiter.shut_down_gnina()

cmd = "gnina -r /home/qzj517/POR-DD/data/raw_data/por_structures/3ES9_1_reduced.pdb -l /home/qzj517/POR-DD/data/raw_data/cyp_screen/first_couple_ligs.sdf --minimize --continuous_operation -o gnina_out.sdf"

asyncio.run(idk(cmd))

